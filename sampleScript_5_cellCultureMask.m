% This script demonstrates how to a) load a stack of data stored in a .mat
% file, rather than in image files, b) apply binary mask to it to discard
% image background. In this case, this is applied to remove the
% surroundings of a Petri dish with a cell culture. The script ultimately
% plots the conduction map (based on calcium mapping data).

clear parameters

data = load('data/cultureStack/stack.mat');
imStackCulture = data.imStackCulture;

% This is what the average frame looks like - it is stored and then its
% copy can be painted over with black color in the areas that are not to be
% processed.
imForMaskCulture = COSMAS.getImageForMask(uint16(imStackCulture));
imwrite(imForMaskCulture, 'masks/cultureAvg.png');

% Now we applied red color to parts of the image we want to ignore
% (assuming the rest is grayscale) - binary mask is then obtained by
% keeping only pixels that were not painted in red color.
% In this case, the mask was drawn manually, but in general, it can be
% often also constructed automatically based on imForMaskCulture (e.g. via
% intensity thresholding)
mask = imread('masks/cultureMask.png');
maskBinary = mask(:,:,1) == mask(:,:,2); % if red==green (i.e., no red painting), we keep the image.

imStackCultureMasked = COSMAS.applyMask(imStackCulture, maskBinary);

parameters.spikesPointDown = false;
parameters.baselineDefinition = 'first';
parameters.smoothingParameter = 21;
[baseline, amplitude, duration, activationMaps, clockFiltered] = COSMAS.analyseRegularPacing(imStackCultureMasked, 200, parameters); % 1 Hz, but 200 fps recording, hence bcl = 200.
am =activationMaps.mapMean;
figure(15);
contourf(flipud(am - quantile(am(:), 0.01)), 'levels', 1, 'ShowText', 'on'); caxis([0 12]);
