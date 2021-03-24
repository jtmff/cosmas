% This script takes two recordings of voltage mapping in the same heart,
% before and after application of isoproterenol, showing an increase in
% conduction velocity under isoproterenol.
clear parameters params
% Reading png images in the 'controlV' folder, performing 2x2 binning
paramReading.extension = '.png';
paramReading.binningFactor = 4;
imStackV = COSMAS.readFolder('data/controlV', paramReading);

% Setting parameters of processing
parameters.spikesPointDown = true;
parameters.baselineDefinition = 'firstlast';
parameters.objectDetection = 'augmented';
bcl = 150; % a stimulus was applied every 150 ms

% Analysing the data
[base, ampl, dur, activationMaps, rm, cf] = COSMAS.analyseRegularPacing(imStackV, bcl, parameters);

% Plotting the mean activation map and local CV
COSMAS.plotActivationMap(activationMaps.mapMean);
COSMAS.getLocalCV(activationMaps.mapMean, 4, 5, 2); 

% reading and analysing the folder with recording of the heart with
% isoproterenol
[imStackViso, at] = COSMAS.readFolder('data/isoV', paramReading);
[baselineISO, amplitudeISO, durationISO, activationMapsISO, recoveryMapsISO, clockFilteredISO] = COSMAS.analyseRegularPacing(imStackViso, bcl, parameters);

% taking mean activation maps for control and iso
am = activationMaps.mapMean;
amISO = activationMapsISO.mapMean;

figure(3);
subplot(1,2,1)
contourf(flipud(am - min(min(am))), 'levels', 1, 'ShowText', 'on');
title('Control');
subplot(1,2,2)
contourf(flipud(amISO - min(min(amISO))), 'levels', 1, 'ShowText', 'on');
title('ISO');

% Measuring conduction velocity between row=13,column=2 and row=1,column=2.
% The factor of 7/16 is the pixel side length in mm (FOV is 7-by-7mm and we
% have 16-by-16 pixels after the binning), and 1000 is the sampling
% frequency (1000 Hz).
cv1 = COSMAS.getCV(am, [13, 2, 1, 2], [7/16, 1000]); 
cv2 = COSMAS.getCV(amISO, [13, 2, 1, 2], [7/16, 1000]);
disp(['Relative speedup: ' num2str(cv2/cv1)]);
