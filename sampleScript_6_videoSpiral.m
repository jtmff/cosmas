% This scripts shows how to a) read data of a spiral wave mapping using
% off-axis illumination imaging from an mp4 video, b) plot the activation
% map and show how standard deviation of duration may be used to get rid of
% spurious signal areas.

clear 

binningFactor = 4; % we do 4x4 binning on the unnecessarily high-res video.

% reading the video - takes a while.
v = VideoReader('data/spiral_wave.mp4');
stack = zeros(v.Height/binningFactor, v.Width/binningFactor, v.NumFrames);
for i = 1:v.NumFrames
%     disp(i);
    frame = v.readFrame();
    frame = rgb2gray(frame);
    frame = COSMAS.binningVec(frame, log2(binningFactor));
    stack(:,:,i) = frame;
end

%%
parameters.baselineSubtractionOrder = 3;
parameters.smoothingParameter = 21;
parameters.spikesPointDown = false;
parameters.baselineDefinition = 'firstlast';
parameters.objectDetection = 'first';

bcl = 34; % in frames 
msPerFrame = 10; % culture imaged at 100 fps

% Feature extraction is performed (given noisy signal and pixels with
% irregular activaton, there are some instances of 'issue in getting half
% activation' - this is just a warning, not a problematic error. The
% resulting pixels will have the feature values 'NaN'.
[baseline, amplitude, duration, activationMaps, recoveryMaps, clockV] = COSMAS.analyseRegularPacing(stack, bcl, parameters);

% activation map of the spiral wave is plotted
COSMAS.plotActivationMap(activationMaps.mapMean * msPerFrame, 1, 20); % a contour line is drawn every 20 ms
set(gcf,'Position',[100 100 800 800]);

%% Making a mask for visualisation cleaning based on postprocessing of duration maps
% (a manually drawn mask might be more accurate nevertheless, this is only
% a proof of concept)

% The standard deviation is computed over the 3rd dimension of the stack of
% signal durations, showing high values in areas of irregular activation,
% while giving low values where the main spiral is present.
mapStd = nanstd(duration.maps*msPerFrame, 0, 3);

figure(2); clf;
subplot(1,3,1);
imagesc(mapStd); colorbar; title('Std of duration (over wave passes) map');

maskToDiscard = isnan(mapStd) | (mapStd > 3*msPerFrame); % a fairly arbitrary threshold

subplot(1,3,2); 
imshow(maskToDiscard); title('for removal - before filtering');

% morphological operations
str = strel('sphere', 1);
strLarger = strel('sphere', 3);
maskToDiscard = imfill(maskToDiscard,'holes');
maskToDiscard = imopen(maskToDiscard, str);
maskToDiscard = imclose(maskToDiscard, strLarger);

subplot(1,3,3);
imshow(maskToDiscard); title('for removal - after filtering');
set(gcf,'Position',[100 100 1800 500])

% filtering the activation map to set 'for removal' pixels to NaN
actMapFiltered = activationMaps.mapMean;
actMapFiltered(maskToDiscard) = NaN;

% plotting activation map after postprocessing
COSMAS.plotActivationMap(actMapFiltered*msPerFrame, 3, 20);