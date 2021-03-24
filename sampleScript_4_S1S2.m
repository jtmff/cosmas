% This script demonstrates how to use custom comb to process S1S2
% recordings (3 S1 stimuli, 1 S2 stimulus), producing a calcium transient restitution curve (Figure 6 in
% the article). For the purposes of this script, we will take the average
% of the field of view as a single trace on which we measure the
% restitution.
% In addition, it is demonstrated at the end of the code how the processing may be ran using
% parallel for-loop, accelerating the runtime.
clear parameters params
tic
s2s = [60:10:130]; % the range of s2 intervals we will process (we know that s1=150 ms)

for iS2 = 1:length(s2s)
    % Read the data using 2x2 binning
    paramReading.extension = '.png';
    paramReading.binningFactor = 2;
    [imStackCa, avgTrace] = COSMAS.readFolder(['data/caRestitution/p' num2str(s2s(iS2)) '.ma/'], paramReading);
    
    
    stackAvgTrace = COSMAS.traceToStack(avgTrace); % we convert the average trace of the recording (a vector of 1000 x 1 elements) to a "fake" 3D stack (1 x 1 x 1000), so that it can be seamlessly used in COSMAS.analyseRegularPacing
    stackAvgTrace = stackAvgTrace(1:700); % we remove the last 300 ms as there is not anything interesting there.
    sat{iS2} = avgTrace; % storing the trace for future plotting
    
    % Setting parameters of processing and processing itself
    parameters.spikesPointDown = false;
    parameters.baselineDefinition = 'first';
    parameters.objectDetection = 'largest';
    
    % A custom comb is defined, so that:
    % a) the first diastole (Ca signal minimum) is to be found between 1 and 30th frame (as
    % given by the first element)
    % b) the second and third diastole come after further 150 and then
    % another 150 ms (because we know that S1 = 150 ms)
    % c) the next diastole (after the 3rd S1, so between the last S1 and
    % the S2) comes after the s2s(iS2) ms (60,70,80,..., depending on where
    % we are in the for loop)
    % d) And the next diastole (after the S2) comes after 100 more frames.
    % This is not critical in general, but in case there are some strange
    % phenomena after the S2, using this "extra" diastole means, that we
    % get quite tightly segmented S2 stimulus and we can just discard the
    % last processed interval (from the diastole after the S2 to the end of
    % the recording).
    %
    % - cumsum is used so that the comb can be interpreted as indices of
    % the comb tips that are then shifted over the signal. In this way, the
    % first element furthermore encodes up to which frame we can slide such
    % a comb. 
    parameters.customComb = cumsum([30,150,150,s2s(iS2),100]);  

    bcl = 150; % this is not important in this script, as the custom comb is used

    [baseline, amplitude, duration, activationMaps, recoveryMaps, clockFiltered] = COSMAS.analyseRegularPacing(stackAvgTrace, bcl, parameters);
    
    relativeCaT(iS2) = amplitude.data(4)/mean(amplitude.data(1:3)); % to get relative CaT amplitude during S2 stimulus, we divide the absolute amplitude by the average of the CaT amplitude of the three S1 stimuli.
    %- note that amplitude.data contains one more element (NaN) - that
    %corresponds to what is described above in the point (d), which is the
    %ignored part of the recording.
end


%% Plotting

figure(1); clf;
plot(s2s, relativeCaT,'d-' , 'LineWidth', 2);
xlabel('S2 (ms)', 'FontSize',14);
ylabel('Relative Ca release', 'FontSize',14);
set(findall(gcf,'-property','FontSize'),'FontSize',14)

figure(2); clf
plot(sat{1}(1:700), 'LineWidth', 1.5);
xlabel('Time (ms)', 'FontSize',14);
ylabel('Signal intensity', 'FontSize',14);
set(findall(gcf,'-property','FontSize'),'FontSize',14)

figure(3); clf
plot(sat{8}(1:700), 'LineWidth', 1.5);
xlabel('Time (ms)', 'FontSize',14);
ylabel('Signal intensity', 'FontSize',14);
set(findall(gcf,'-property','FontSize'),'FontSize',14)


%% How to run COSMAS in parallel - uncomment the code below and replace the first code cell with this. It is similar to the simple for-loop above, but because
% of how Matlab handles parallelism, paramReading and parameters are now
% arrays, with each element corresponding to one s2 protocol

% clear parameters params
% tic
% s2s = [60:10:130];
% parfor iS2 = 1:length(s2s)
%     paramReading(iS2).extension = '.png';
%     paramReading(iS2).binningFactor = 2;
%     [imStackCa, avgTrace] = COSMAS.readFolder(['caRestitution/p' num2str(s2s(iS2)) '.ma/'], paramReading(iS2));
%     stackAvgTrace = COSMAS.traceToStack(avgTrace);
%     stackAvgTrace = stackAvgTrace(1:700); % we remove the empty end of the recording.
%     sat{iS2} = avgTrace;
%     % Setting parameters of processing and processing itself
%     parameters(iS2).spikesPointDown = false;
%     parameters(iS2).baselineDefinition = 'first';
%     parameters(iS2).objectDetection = 'largest';
%     parameters(iS2).customComb = cumsum([30,150,150,s2s(iS2),100]); 
% 
%     bcl = 150; % unimportant, because custom comb is used
% 
%     [baseline, amplitude, duration, activationMaps, recoveryMaps, clockFiltered] = COSMAS.analyseRegularPacing(stackAvgTrace, bcl, parameters(iS2));
%     
%     relativeCaT(iS2) = amplitude.data(4)/mean(amplitude.data(1:3));
% end
% 
