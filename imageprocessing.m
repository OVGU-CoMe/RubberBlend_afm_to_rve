function [filler_content_A,filler_content_B,binary_image] = imageprocessing(start_image,blend_ratio,filler_content,interphase_width)

%% Description
% author: IWTM @ OVGU
% date:   2024/11/18
%
% This function processes AFM images of filled binary elastomer blends,
% e.g. carbon black-filled NR-SBR. It prepares a binary image as an RVE
% (representative volume element) for FEA (finite element analysis) and
% determines the filler content of the polymer phases. The interphase
% between the polymer phases is obtained from blurring a sharp morphology.
%
% input: % start_image (rgb, cut to size)
%        % blend_ratio (fraction of the "white" phase, here NR)
%        % filler_content (overall filler content in phr)
%        % interphase_width
%
% output: % filler_content_A (filler content of the "white" phase, here NR)
%         % filler_content_B (filler content of the darker phase, here SBR)
%         % binary_image (white = phase A, black = phase B)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare parameters

% absolute error for determining threshold values
fehler = 0;

% polymer fractions
fraction_white = blend_ratio;
fraction_gray  = 1- blend_ratio;

% filler content (volume-related quantity)
filler_volume = filler_content/(200+filler_content);

%% Create Three-Phases Image (phase A, phase B, filler)

% transform image to gray
Igray = rgb2gray(start_image);

% determine threshold values between gray and black (phase B and filler)
% and define mask for black (filler)
thresholdvalue1 = 1;
number_total_pixels = numel(Igray);
differenz = 1;
while differenz > fehler
    mask_black = Igray < thresholdvalue1;
    number_black_pixels = nnz(mask_black);
    area_black = number_black_pixels/number_total_pixels;   % filler
    differenz = filler_volume - area_black;
    thresholdvalue1 = thresholdvalue1 + 1;
end
mask_black = Igray < thresholdvalue1;
number_black_pixels = nnz(mask_black);
area_black = number_black_pixels/number_total_pixels;   % filler
differenz_neu = filler_volume - area_black;
% check whether previous are current step is better
if abs(differenz) < abs(differenz_neu)
    thresholdvalue1 = thresholdvalue1 -1;
    mask_black = Igray < thresholdvalue1;
    number_black_pixels = nnz(mask_black);
    area_black = number_black_pixels/number_total_pixels;   % filler
end

% determine threshold values between gray and white (phase B and phase A)
% and define mask for gray (phase B)
thresholdvalue2 = 255;
differenz       = 1;
while differenz > fehler
    mask_gray = Igray > thresholdvalue1 & Igray <= thresholdvalue2;
    number_gray_pixels  = nnz(mask_gray);
    area_gray  = number_gray_pixels/number_total_pixels;    % SBR
    area_white = 1- area_gray - area_black;                 % NR
    fraction_p_gray = area_gray/(area_gray + area_white);

    differenz = fraction_p_gray - fraction_gray;
    thresholdvalue2 = thresholdvalue2 - 1;
end
mask_gray = Igray > thresholdvalue1 & Igray <= thresholdvalue2;
number_gray_pixels  = nnz(mask_gray);
area_gray  = number_gray_pixels/number_total_pixels;    % SBR
area_white = 1- area_gray - area_black;                 % NR
fraction_p_gray = area_gray/(area_gray + area_white);
differenz_neu = fraction_p_gray - fraction_gray;
% check whether previous are current step is better
if abs(differenz) < abs(differenz_neu)
    thresholdvalue2 = thresholdvalue2 + 1;
    mask_gray = Igray > thresholdvalue1 & Igray <= thresholdvalue2;
    number_gray_pixels  = nnz(mask_gray);
    area_gray  = number_gray_pixels/number_total_pixels;    % SBR
    area_white = 1- area_gray - area_black;                 % NR
    fraction_p_gray = area_gray/(area_gray + area_white);
end

% define mask for white (phase A)
mask_white = Igray > thresholdvalue2;
fraction_p_white = area_white/(area_gray + area_white);

% initialize resulting image
resultImage = uint8(zeros(size(Igray)));

% set pixels
resultImage(mask_white) = 255;
resultImage(mask_gray) = 128;
resultImage(mask_black) = 0;

% save three-phases image as file
imwrite(resultImage,'three_phases.png')

%% Determine filler contents in both polymer phases

% obtain boundaries of filler particles
boundaries_black = bwperim(mask_black,8);
SE = strel('disk',5);
boundaries_black = imdilate(boundaries_black,SE);

% define overlap between boundaries and white or gray masks
overlap_white = boundaries_black & mask_white;
overlap_gray = boundaries_black & mask_gray;

% count pixels
pixels_overlap_white = nnz(overlap_white);
pixels_overlap_gray  = nnz(overlap_gray);
pixels_overlap = pixels_overlap_white + pixels_overlap_gray;

% compute fractions from the pixels
fraction_surrounding_white = pixels_overlap_white / pixels_overlap;
fraction_surrounding_gray  = pixels_overlap_gray / pixels_overlap;

% convert fraction to filler content in phr in both phases
filler_content_white = filler_content * fraction_surrounding_white / fraction_white;
filler_content_gray = filler_content * fraction_surrounding_gray / fraction_gray;

%% Create binary image (without filler phase)

% find out which phase to grow
if fraction_surrounding_white < fraction_surrounding_gray
    overlap_to_grow = overlap_white;
    fraction_to_grow = fraction_surrounding_white;
    value_grown = 255;
    other_value = 128;
else
    overlap_to_grow = overlap_gray;
    fraction_to_grow = fraction_surrounding_gray;
    value_grown = 128;
    other_value = 255;
end

% grow the overlap into the gaps
% overlap wachsen lassen
fraction_grown = 0.0;
SE = strel('disk',2);
% loop for successive growing
while fraction_grown < fraction_to_grow
    mask_grown = imdilate(overlap_to_grow,SE); % grow
    mask_grown_in_black = mask_grown & mask_black; % how much has grown?
    fraction_grown = nnz(mask_grown_in_black)/ number_black_pixels;
    overlap_to_grow = mask_grown_in_black; % prepare next iteration
end

% construct two-phases image
result_two_phases = uint8(zeros(size(Igray)));
result_two_phases(mask_white) = 255;
result_two_phases(mask_gray) = 128;
result_two_phases(mask_grown_in_black) = value_grown;
remaining_black =~(mask_white|mask_gray|mask_grown_in_black);
result_two_phases(remaining_black) = other_value;

% save two-phases image as file
imwrite(imbinarize(result_two_phases),'two_phases.png')
%% Create interphase
% keep two phases image if width is zero
if interphase_width == 0
    Iblur = imbinarize(result_two_phases);

% blur two-phases image 
else
    Iblur = imgaussfilt(double(imbinarize(result_two_phases)),interphase_width);
end

%% write output
% filler contents
filler_content_A = filler_content_white;
filler_content_B = filler_content_gray;

% binary image
imwrite(Iblur,'binary_image.png')
binary_image =Iblur;