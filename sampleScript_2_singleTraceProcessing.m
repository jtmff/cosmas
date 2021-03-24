% This script shows how to extract properties such as APD from single-trace
% recordings. While such a trace here was extracted as a spatial average of
% imaging experiment, electrophysiological recordings may be naturally
% processed in a similar fashion.
clear parameters params
paramReading.extension = '.png';
paramReading.binningFactor = 2;
% the second output parameter of readFolder is the signal trace produced as a
% result of spatially averaging the field of view in each frame.
[imStackViso, at] = COSMAS.readFolder('data/isoV', paramReading);

% the trace (a vector) is converted to a "3D stack" with dimensions
% 1-by-1-by-n so that analyseRegularPacing may be directly applied
stackTrace = COSMAS.traceToStack(at);
parameters.spikesPointDown = true; % analysing membrane potential, with dye intensity reducing upon membrane depolarization
parameters.baselineDefinition = 'firstlast';

% analysing the trace
[baseline, amplitude, duration, activationMaps, clockFiltered] = COSMAS.analyseRegularPacing(stackTrace, 150, parameters);

% plotting the trace and showing the mean APD in the title
% (Please mind that taking average APD in a spatially averaged recording
% has only limited relevance and it may be often more representative to
% measure APD for each pixel and then average that - especially for preps
% with slow conduction such as certain cell cultures, the APD of
% whole-FOV-average is strongly modulated by conduction velocity. While a
% recording with infinite conduction velocity would have the same
% single-pixel and whole-FOV APD, if the conduction velocity is low, the
% temporally shifted action potentials combine into one that is possibly
% substantially prolonged compared to the single-pixel one).
figure(1); clf;
plot(at); 
xlabel('Frame (no.)'); 
ylabel('voltage dye intensity');
title(['APD75 observed in the recording: ' num2str(duration.dataMean)]); % APD75 is the default level for APD extraction.
