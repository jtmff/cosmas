
% TODO copyright stuff
% TODO copy/mention licence of natsort and natsortfiles
classdef COSMAS
    % COSMAS A matlab library for processing of optical mapping data
    
    methods (Static )
        %%%%% READING DATA %%%%%
        function [imStack, avgTrace] = readFolder(folderName, varargin)
            % a helper function returning a stack of images from folder with
            % image sequence. Uses natsort to get natural ordering (this is useful when there are
            % files 1.png, 2.png,...,10.png,11.png, which would get sorted as 1,10,11,2, in purely alphabetic sorting, which is not the right order of frames)
            %
            % IN:
            % folderName is the folder containing a sequence of images
            % corresponding to frames.
            %
            % varargin can contain a structure of parameters (called e.g. readerParameters), with the
            % following fields:
            % 1) readerParameters.extension - extension of images in the
            % given folder.
            % 2) readerParameters.fromTo - a 2-by-1 vector [from to],
            % determining from which to which frame is the recording to be
            % processed.
            % 3) readerParameters.binningFactor - spatial binning parameter
            % (must be a power of two) - binningFactor x binningFactor
            % pixels are aggregated into a single one.
            %
            % OUT:
            % imStack - the resulting 3D stack, where rows and columns
            % correspond to frame rows and columns, and 3rd dimension to
            % index of image in the sequence (i.e. time).
            %
            % avgTrace - a trace corresponding to spatially averaged stack.
            
            
            
            %% Setting defaults, reading parameters, getting filenames to be read.
            extension = '.tif';
            fromTo = [];
            binningFactor = 1;
            if ~isempty(varargin)
                readerParameters = varargin{1};
            end
            
            if (isfield(readerParameters, 'extension'))
                extension = readerParameters.extension;
            end
            
            if (isfield(readerParameters, 'fromTo'))
                fromTo = readerParameters.fromTo;
            end
            
            if (isfield(readerParameters, 'binningFactor'))
                binningFactor = readerParameters.binningFactor;
            end
            
            fnames = dir([folderName '/*' extension]);
            fnames = struct2cell(fnames);
            fnames = fnames(1,:);
            fnames = COSMAS.natsortfiles(fnames);
            
            if isempty(fromTo)
                fromTo = [1 length(fnames)];
            end
            
            %% Checking parameter correctness
            assert(fromTo(1)>=1 & fromTo(2) <= length(fnames), 'Boundaries from-to the frame are incorrect. Fewer images are available than the boundaries ask for.')
            assert(round(log2(binningFactor)) == log2(binningFactor), 'The binning factor must be a power of 2.')
            
            %% Reading stack
            im = imread([folderName '/' fnames{1}]);
            nFrames = fromTo(2) - fromTo(1) + 1;
            
            imStack = zeros(size(im,1)/binningFactor, size(im,2)/binningFactor, nFrames);
            avgTrace = zeros(nFrames, 1);
            
            iFrom = fromTo(1);
            iTo = fromTo(2);
            
            for iFrame = iFrom:iTo
                img = imread([folderName '/' fnames{iFrame}]);
                if (size(img,3) == 3)
                    img = rgb2gray(img);
                end
                
                if (binningFactor ~= 1)
                    img = COSMAS.binningVec(img, log2(binningFactor));
                end
                imStack(:,:,iFrame - iFrom + 1) = img;
                avgTrace(iFrame - iFrom + 1) = nanmean(img(:));
            end
        end
        
        function [imStack, varargout] = readTifStack(fname, varargin)
            % returning a stack of images from a single tif stack.
            % fname is the path to the file (including filename and extension).
            %
            % IN:
            % varargin can contain a structure of parameters (called readerParameters), with the
            % following fields:
            % 1) readerParameters.fromTo - a 2-by-1 vector [from to],
            % determining from which to which frame is the recording to be
            % processed
            % 2) readerParameters.binningFactor - spatial binning parameter
            % (must be a power of two) - binningFactor x binningFactor
            % pixels are aggregated into a single one.
            %
            % OUT:
            % imStack - the resulting 3D stack.
            %
            % avgTrace - a trace corresponding to spatially averaged stack.
            
            % often there is an annotation missing in TIFF stacks, giving
            % annoying warnings, so we turn that off
            warning('off', 'imageio:tifftagsread:noTypeFormatId');
            warning('off', 'imageio:tifftagsread:expectedTagDataFormat');
            info = imfinfo(fname);
            warning('on', 'imageio:tifftagsread:noTypeFormatId');
            warning('on', 'imageio:tifftagsread:expectedTagDataFormat');
            
            nFramesTotal = numel(info);
            
            fromTo = [];
            binningFactor = 1;
            if ~isempty(varargin)
                readerParameters = varargin{1};
            end
            
            if (isfield(readerParameters, 'fromTo'))
                fromTo = readerParameters.fromTo;
            end
            
            if (isfield(readerParameters, 'binningFactor'))
                binningFactor = readerParameters.binningFactor;
            end
            
            if isempty(fromTo)
                fromTo = [1 nFramesTotal];
            end
            
            %% Checking parameter correctness
            assert(fromTo(1)>=1 & fromTo(2) <= nFramesTotal, 'Boundaries from-to the frame are incorrect. Fewer images are available than the boundaries ask for.')
            assert(round(log2(binningFactor)) == log2(binningFactor), 'The binning factor must be a power of 2.')
            
            
            %% Reading stack
            nFrames = fromTo(2) - fromTo(1) + 1;
            imStack = zeros(info(1).Height/binningFactor, info(1).Width/binningFactor, nFrames);
            avgTrace = zeros(nFrames, 1);
            
            iFrom = fromTo(1);
            iTo = fromTo(2);
            
            for iFrame = iFrom:iTo
                img = imread(fname, iFrame, 'Info', info);
                if (size(img,3) == 3)
                    img = rgb2gray(img);
                end
                if (binningFactor ~= 1)
                    img = COSMAS.binningVec(img, log2(binningFactor));
                end
                imStack(:,:,iFrame - iFrom + 1) = img;
                avgTrace(iFrame - iFrom + 1) = nanmean(img(:));
            end
            varargout{1} = avgTrace;
            %% NOTE: the way of reading Tif stacks above works fine and allows reading just the images between selected indices, but it can be quite slow for very large
            % stacks (e.g. 2k frames) - in that case, a faster way (which
            % nevertheless requires going through all the images until iTo)
            % is the one below, available in non-ancient versions of
            % Matlab. One issue there is that if the stack has more than
            % 2^16 frames, this crashes because it uses 16bit counter that
            % overflows :/
            %             tstack = Tiff([folderIn '/' fnameStack.name]);
            %             im = tstack.read();
            %             imStack = zeros(size(im, 1), size(im, 2),length(stackInfo));
            %             for i = 2:length(stackInfo)
            %                 img = tstack.nextDirectory();
            %                 if (binningFactor ~= 1)
            %                     img = COSMAS.binningVec(img, log2(binningFactor));
            %                 end
            %                 imStack(:,:,i) = tstack.read();
            %             end
        end
        
        
        %%%%% ANALYSING DATA %%%%%
        function [baseline, amplitude, duration, activationMaps, recoveryMaps, recordingClock] = analyseRegularPacing(imageStack, bcl, parameters)
            % A function for processing stacks corresponding to
            % recordings with multiple passes of a wave.
            %
            % IN
            % imageStack - a stack representing the recording.
            %
            % bcl - basic cycle length of the recording (in frames).
            %
            % parameters - a structure with parameters:
            %
            % parameters.baselineDefinition - whether baseline of an
            % activation (e.g. calcium transient) is taken as the first
            % element in each segmented activation ('first'), or as the
            % average of the first and last element ('firstlast').
            %
            % parameters.baselineSubtractionOrder - order of polynomial
            % subtraction of signal baseline. Use -1 if no baseline
            % subtraction is to be done (or just leave the parameter field
            % undefined).
            %
            % parameters.durationLevel - level at which duration is
            % extracted, scaled to 0-1. E.g., use 0.8 for APD80.
            %
            % parameters.objectDetection - when multiple activations are discovered in
            % a single signal segment (while only one can be true calcium
            % transient/action potential), this parameter determines how
            % the correct activation is detected (hopefully :)). 'first' -
            % first object found in the segment. 'largest' - the largest
            % object is picked (with most frames). 'augmented' -
            % information on derivative of the signal is used, see the
            % publication. For voltage mapping with multiple wave passes,
            % we recommend using 'augmented' (which is not great for calcium,
            % given that there are no sharp upstrokes). Otherwise, 'largest'
            % is usually more robust than 'first'.
            %
            % parameters.smoothingParameter - width of Savitzky-Golay
            % filtering for signal smoothing. It should be an odd number
            % (if an even number is given, 1 will be added). Given that 4th
            % order smoothing is used, this parameter has to be at least 5
            % when provided (when smaller, no filtering is done).
            %
            % parameters.spikesPointDown - if true, signal activation
            % manifests as reduction in signal intensity (e.g. some voltage
            % dyes), i.e., action potentials "point down". default =
            % false.
            %
            % parameters.verbose - if true, the code reports when there is
            % a problem with segmentation and/or processing of a pixel
            % trace (which can happen when a pixel contains only noise, for
            % example).
            %
            % parameters.waveProcessing – if 'perbeat', each wave is processed separately,
            % and the output structures contain a map for each complete wave pass.
            % If 'hybrid', a recording clock is used to chop the recording into sub-stacks, w
            % hich are then averaged, and a single-wave processing is applied to this subsequently.
            % The value 'hybrid' is good for very noisy recordings and activation mapping,
            % but is not suggested to be used for APD mapping or amplitude measurements.
            % If ‘hybrid’ is used, parameter.objectDetection should be set to ‘largest’.
            % The default is ‘perbeat’.
            %
            % OUT:
            % baseline - a structure describing signal baseline (e.g. bases
            % of calcium transients
            %
            % amplitude - a structure describing signal amplitude
            %
            % duration - a structure describing signal duration (e.g. APD)
            %
            % activationMaps - a structure describing activation pattern in
            % the recording (relative to recording clock; if you want to
            % have minimum activation in 0, just subtract minimum of these maps)
            %
            % recoveryMaps - a structure describing recovery pattern (e.g.
            % the time when APD80 is reached). Relative to the same clock
            % as activationMaps
            %
            % recordingClock - the recording clock determining global
            % synchronization of single segmented activations
            %
            % the structures baseline, amplitude, and duration have the
            % following fields (shown for duration):
            %    duration.maps - a 3D stack where each slice corresponds to
            %    a map of the feature in a single wave pass
            %
            %    duration.mapMean - the average of the previous maps
            %
            %    duration.data - a vector of spatial averages of maps in duration.maps
            %
            %    duration.dataMean = nanmean(duration.data) - the mean of
            %    the field data (i.e., this gives one number summarizing
            %    the whole recording)
            %
            %    duration.mapMeanEven - mean map for even beats (useful for
            %    inspection of alternans).
            %
            %    duration.mapMeanOdd - mean map for odd beats
            %
            %    duration.dataMeanEven - spatial averages of even maps of
            %    the feature
            %
            %    duration.dataMeanOdd - spatial averages of odd maps of the
            %    feature
            %
            % the structure activationMaps and recoveryMaps have the same
            % fields, except the ones starting with 'data' (there is not
            % much point in spatially averaged activation).
            
            %% Reading and verifying parameters
            spikesPointDown = false; % if spikes point down instead of up. This is relevant for finding diastole, and for extracting baseline of the signal
            if (isfield(parameters, 'spikesPointDown'))
                spikesPointDown = parameters.spikesPointDown;
            end
            
            baselineDefinition = 'first';
            if (isfield(parameters, 'baselineDefinition'))
                baselineDefinition = parameters.baselineDefinition;
            end
            assert(ismember(baselineDefinition,{'first','firstlast'}), 'baselineDefinition must be either first or firstlast');
            
            objectDetection = 'largest';
            if (isfield(parameters, 'objectDetection'))
                objectDetection = parameters.objectDetection;
            end
            assert(ismember(objectDetection,{'augmented','first','largest'}), 'objectDetection must be either augmented, first or largest');
            
            
            baselineSubtractionOrder = -1;
            if (isfield(parameters, 'baselineSubtractionOrder'))
                baselineSubtractionOrder = parameters.baselineSubtractionOrder;
            end
            
            smoothingParameter = 11;
            if (isfield(parameters, 'smoothingParameter'))
                smoothingParameter = parameters.smoothingParameter;
            end
            
            durationLevel = 0.75;
            if (isfield(parameters, 'durationLevel'))
                durationLevel = parameters.durationLevel;
            end
            
            verbose = false;
            if (isfield(parameters, 'verbose'))
                verbose = parameters.verbose;
            end
            assert(islogical(verbose), 'verbose parameter must be either true or false');
            
            waveProcessing = 'perbeat';
            if (isfield(parameters, 'waveProcessing'))
                waveProcessing = parameters.waveProcessing;
            end
            assert(ismember(waveProcessing,{'perbeat','hybrid'}), 'waveProcessing must be perbeat or hybrid');
            
            customComb = [];
            if (isfield(parameters, 'customComb'))
                customComb = parameters.customComb;
            end
            
            nRows = size(imageStack, 1);
            nCols = size(imageStack, 2);
            %% Getting global signal, from which we extract the recording "timer" that separates action potentials/calcium transients
            % Unlike usual, we extract locations of peaks p1,p2,..., and then build a
            % timer that starts at (p1+p2)/2, so peaks happen roughly at
            % midpoint between the points in the timer
            avgTrace = squeeze(mean(mean(imageStack,1),2));
            
            if (baselineSubtractionOrder > 0)
                [~, avgTrace] = COSMAS.smoothTrace(avgTrace, smoothingParameter, baselineSubtractionOrder);
            else
                avgTrace = COSMAS.smoothTrace(avgTrace, smoothingParameter);
            end
            
            if (spikesPointDown)
                maxActivations = COSMAS.combGetMinima(avgTrace, bcl, [], customComb); % as a 2nd parameter, we pass a pair of empty vector (we don't vary the surroundings of comb tips here), and the custom comb, if provided.
            else
                maxActivations = COSMAS.combGetMinima(-avgTrace, bcl, [], customComb);
            end
            
            recordingClock = mean([maxActivations(1:end-1);maxActivations(2:end)]);
            
            if (~isempty(customComb))
                recordingClock = [1 recordingClock size(imageStack, 3)];
            end
            
            if(strcmp(waveProcessing, 'hybrid')) % hybrid processing, where the timer is used to chop the recording to smaller parts and it is then averaged and processed using singleWavePass
                recordingClock = recordingClock(1):bcl:length(avgTrace);
                
                % Now, we process all traces in the stack with smoothing
                % and baseline subtraction, then chopping it to bcl-sized
                % chunks and averaging them, before processing it as a
                % single wave
                for iRow = 1:nRows
                    for iCol = 1:nCols
                        % extracting and smoothing the trace
                        pixelTrace = squeeze(imageStack(iRow,iCol,:));
                        if (sum(isnan(pixelTrace))>0)
                            continue;
                        end
                        
                        if (baselineSubtractionOrder > 0)
                            [~, traceSmoothed] = COSMAS.smoothTrace(pixelTrace, smoothingParameter, baselineSubtractionOrder);
                        else
                            traceSmoothed = COSMAS.smoothTrace(pixelTrace, smoothingParameter);
                        end
                        imageStack(iRow, iCol, :) = traceSmoothed;
                    end
                end
                
                stackAvg = zeros(nRows, nCols, bcl);
                for iPass = 1:(length(recordingClock) - 1)
                    stackAvg = stackAvg + imageStack(:,:,recordingClock(iPass):(recordingClock(iPass+1)-1));
                    traces{iPass} = avgTrace(recordingClock(iPass):(recordingClock(iPass+1)-1));
                end
                stackAvg = stackAvg./(length(recordingClock) - 1);
                
                
                [baseline, amplitude, duration, activationMaps, recoveryMaps] = COSMAS.analyseSinglePass(stackAvg, parameters);
                return;
            else % otherwise we do standard processing, beat per beat
                
                mapsBaseline = nan(nRows, nCols, length(recordingClock)-1); % for each pixel, signal baseline in the i-th pass of the wave
                mapsAmplitude = nan(nRows, nCols, length(recordingClock)-1); % for each pixel, signal amplitude in the i-th pass of the wave
                mapsDuration = nan(nRows, nCols, length(recordingClock)-1); % for each pixel, duration in the i-th pass of the wave
                mapsActivationTimes = nan(nRows, nCols, length(recordingClock)-1); % and activation times
                mapsRecoveryTimes = nan(nRows, nCols, length(recordingClock)-1); % and activation times
                
                % for each trace, we segment it using comb, and use local
                % maxima to assign the found values using recordingClock. We
                % don't want to use recordingClock or anything like that itself
                % for segmentation of spikes/transients, as that does not give
                % fine-enough information (e.g. in discordant alternans)
                for iRow = 1:nRows
                    for iCol = 1:nCols
                        pixelTrace = squeeze(imageStack(iRow, iCol, :));
                        
                        if (sum(isnan(pixelTrace))>0)
                            continue;
                        end
                        
                        if (baselineSubtractionOrder > 0)
                            [~, traceSmoothed] = COSMAS.smoothTrace(pixelTrace, smoothingParameter, baselineSubtractionOrder);
                        else
                            traceSmoothed = COSMAS.smoothTrace(pixelTrace, smoothingParameter);
                        end
                        
                        % For the trace, we extract properties between all its
                        % diastoles
                        
                        if (spikesPointDown)
                            diastoles = COSMAS.combGetMinima(-traceSmoothed, bcl, [], customComb);
                        else
                            diastoles = COSMAS.combGetMinima(traceSmoothed, bcl, [], customComb);
                        end
                        
                        % we also find minima/maxima of dv/dt that serve
                        % 'augmented' object detection(finding objects nearest peak
                        % diff(signal).
                        diffSignal = COSMAS.smoothTrace(diff(traceSmoothed), smoothingParameter);
                        
                        
                        if (spikesPointDown)
                            peakDiffs = COSMAS.combGetMinima(diffSignal, bcl, [], customComb);
                        else
                            peakDiffs = COSMAS.combGetMinima(-diffSignal, bcl, [], customComb);
                        end
                        
                        for iStart = 1:(length(diastoles)-1)
                            timeStart = diastoles(iStart);
                            timeEnd = diastoles(iStart + 1);
                            
                            iCenter = (timeStart + timeEnd)/2; % converting the peak within single activation transient to the global temporal coordinates.
                            iBin = sum(iCenter > recordingClock); % after how many elements of the recording clock does the location come?
                            
                            if iBin<1 || iBin>=length(recordingClock) % before/after clock starts
                                continue
                            end
                            
                            peakDiffActivation = peakDiffs(peakDiffs>=timeStart & peakDiffs<=timeEnd) - timeStart + 1; % we take the time of peak activation relative to timeStart
                            
                            % single activation is extracted and processed
                            activationTrace = traceSmoothed(timeStart:timeEnd);
                            [saBaseline, saAmplitude, saDuration, saActivation, saRecovery] = COSMAS.processSingleActivation(activationTrace, spikesPointDown, baselineDefinition,durationLevel, objectDetection,recordingClock, timeStart, peakDiffActivation, verbose);
                            
                            
                            mapsBaseline(iRow, iCol, iBin) = saBaseline;
                            mapsAmplitude(iRow, iCol, iBin) = saAmplitude;
                            mapsDuration(iRow, iCol, iBin) = saDuration;
                            mapsActivationTimes(iRow, iCol, iBin) = saActivation;
                            mapsRecoveryTimes(iRow, iCol, iBin) = saRecovery;
                        end
                    end
                end
                
                %% Now we remove slices with empty entries in activation - this refers to some pixels having not enough information (e.g. not enough time at the end to contain a full action potential/CaT), so these slices are discarded.
                nZerosInSlice = squeeze(sum(sum(mapsActivationTimes==0, 1), 2));
                
                mapsBaseline = mapsBaseline(:, :, nZerosInSlice == 0);
                mapsAmplitude = mapsAmplitude(:, :, nZerosInSlice == 0);
                mapsDuration = mapsDuration(:, :, nZerosInSlice == 0);
                mapsActivationTimes = mapsActivationTimes(:, :, nZerosInSlice == 0);
                
                % returning baseline
                baseline.maps = mapsBaseline;
                baseline.mapMean = nanmean(mapsBaseline, 3);
                baseline.data = squeeze(nanmean(nanmean(mapsBaseline, 1), 2));
                baseline.dataMean = nanmean(baseline.data);
                
                if (size(mapsBaseline, 3) > 1)
                    baseline.mapMeanEven = nanmean(mapsBaseline(:,:,2:2:end), 3);
                    baseline.mapMeanOdd = nanmean(mapsBaseline(:,:,1:2:end), 3);
                    baseline.dataMeanEven = squeeze(nanmean(baseline.data(2:2:end)));
                    baseline.dataMeanOdd = squeeze(nanmean(baseline.data(1:2:end)));
                end
                
                % returning amplitude
                amplitude.maps = mapsAmplitude;
                amplitude.mapMean = nanmean(mapsAmplitude, 3);
                amplitude.data = squeeze(nanmean(nanmean(mapsAmplitude, 1), 2));
                amplitude.dataMean = nanmean(amplitude.data);
                
                if (size(mapsAmplitude, 3) > 1)
                    amplitude.mapMeanEven = nanmean(mapsAmplitude(:,:,2:2:end), 3);
                    amplitude.mapMeanOdd = nanmean(mapsAmplitude(:,:,1:2:end), 3);
                    amplitude.dataMeanEven = squeeze(nanmean(amplitude.data(2:2:end)));
                    amplitude.dataMeanOdd = squeeze(nanmean(amplitude.data(1:2:end)));
                end
                
                
                % returning duration
                duration.maps = mapsDuration;
                duration.mapMean = nanmean(mapsDuration, 3);
                duration.data = squeeze(nanmean(nanmean(mapsDuration, 1), 2));
                duration.dataMean = nanmean(duration.data);
                
                if (size(mapsDuration, 3) > 1)
                    duration.mapMeanEven = nanmean(mapsDuration(:,:,2:2:end), 3);
                    duration.mapMeanOdd = nanmean(mapsDuration(:,:,1:2:end), 3);
                    duration.dataMeanEven = squeeze(nanmean(duration.data(2:2:end)));
                    duration.dataMeanOdd = squeeze(nanmean(duration.data(1:2:end)));
                end
                
                % returning data on activation maps
                activationMaps.maps = mapsActivationTimes;
                if (size(mapsActivationTimes, 3) > 1)
                    activationMaps.mapMean = nanmean(mapsActivationTimes, 3);
                    activationMaps.mapMeanEven = nanmean(mapsActivationTimes(:,:,2:2:end), 3);
                    activationMaps.mapMeanOdd = nanmean(mapsActivationTimes(:,:,1:2:end), 3);
                end
                
                % returning data on recovery maps
                recoveryMaps.maps = mapsRecoveryTimes;
                if (size(mapsRecoveryTimes, 3) > 1)
                    recoveryMaps.mapMean = nanmean(mapsRecoveryTimes, 3);
                    recoveryMaps.mapMeanEven = nanmean(mapsRecoveryTimes(:,:,2:2:end), 3);
                    recoveryMaps.mapMeanOdd = nanmean(mapsRecoveryTimes(:,:,1:2:end), 3);
                end
                
                % returning clock separating activations. This is filtered to
                % remove elements corresponding to incompletely filled maps
                % (i.e. usually this would be just the last map, when it's not
                % completely filled because some pixels don't have a
                % long-enough signal there)
                clockFiltered = recordingClock;
                slicesDeleted = find(nZerosInSlice>0);
                if (~isempty(slicesDeleted))
                    clockFiltered(slicesDeleted + 1) = [];
                end
                recordingClock = clockFiltered;
            end
        end
        
        function [baseline, amplitude, duration, activationMaps, recoveryMaps] = analyseSinglePass(imageStack, parameters)
            % A function for processing stacks corresponding to
            % recordings with single pass of a wave. Please see the
            % documentation of COSMAS.analyseRegularPacing for the
            % description of inputs/outputs.
            %
            % In the output structures, the only fields are mapMean and
            % dataMean. This may look illogical (there is just one map per
            % the recording, so why would one carry out averaging over slices,
            % which is essentially identity?), but
            % it's based on our experience that COSMAS users mainly use the
            % outputs of this function in a similar way as they'd use
            % mapMean/dataMean returned by analyseRegularPacing -
            % therefore, calling the outputs the same means there is less
            % code rewriting needed when one switches between
            % analyseRegularPacing and analyseSinglePass.
            
            spikesPointDown = false; % if spikes point down instead of up. This is relevant for finding diastole, and for extracting baseline of the signal
            if (isfield(parameters, 'spikesPointDown'))
                spikesPointDown = parameters.spikesPointDown;
            end
            
            baselineDefinition = 'first';
            if (isfield(parameters, 'baselineDefinition'))
                baselineDefinition = parameters.baselineDefinition;
            end
            assert(ismember(baselineDefinition,{'first','firstlast'}), 'baselineDefinition must be either first or firstlast');
            
            objectDetection = 'largest';
            if (isfield(parameters, 'objectDetection'))
                objectDetection = parameters.objectDetection;
            end
            assert(ismember(objectDetection,{'augmented','first','largest'}), 'objectDetection must be either augmented, first or largest');
            
            
            baselineSubtractionOrder = -1;
            if (isfield(parameters, 'baselineSubtractionOrder'))
                baselineSubtractionOrder = parameters.baselineSubtractionOrder;
            end
            
            smoothingParameter = 11;
            if (isfield(parameters, 'smoothingParameter'))
                smoothingParameter = parameters.smoothingParameter;
            end
            
            durationLevel = 0.75;
            if (isfield(parameters, 'durationLevel'))
                durationLevel = parameters.durationLevel;
            end
            
            verbose = false;
            if (isfield(parameters, 'verbose'))
                verbose = parameters.verbose;
            end
            assert(islogical(verbose), 'verbose parameter must be either true or false');
            
            nRows = size(imageStack, 1);
            nCols = size(imageStack, 2);
            mapBaseline = nan(nRows, nCols);
            mapAmplitude = nan(nRows, nCols);
            mapDuration = nan(nRows, nCols);
            mapActivationTimes = nan(nRows, nCols);
            mapRecoveryTimes = nan(nRows, nCols);
            
            for iRow = 1:nRows
                for iCol = 1:nCols
                    pixelTrace = squeeze(imageStack(iRow, iCol, :));
                    
                    % we smooth the trace, but do not perform baseline
                    % subtraction - that can cause a huge mess in
                    % single-wave traces
                    traceSmoothed = COSMAS.smoothTrace(pixelTrace, smoothingParameter);
                    
                    
                    if (sum(isnan(pixelTrace))>0)
                        continue;
                    end
                    
                    timeStart = 1;
                    
                    diffSignal = diff(traceSmoothed) + 1;
                    if (spikesPointDown)
                        [~ ,peakDiffActivation] = min(diffSignal);
                    else
                        [~, peakDiffActivation] = min(-diffSignal);
                    end
                    recordingClock = [0, inf]; % this is a fairly dummy value, making sure that the results are considered to belong to the first wave (out of 1 :)).
                    % single activation is extracted and processed
                    [saBaseline, saAmplitude, saDuration, saActivation, saRecovery] = COSMAS.processSingleActivation(traceSmoothed, spikesPointDown, baselineDefinition,durationLevel, objectDetection,recordingClock, timeStart, peakDiffActivation, verbose);
                    
                    mapBaseline(iRow, iCol) = saBaseline;
                    mapAmplitude(iRow, iCol) = saAmplitude;
                    mapDuration(iRow, iCol) = saDuration;
                    mapActivationTimes(iRow, iCol) = saActivation;
                    mapRecoveryTimes(iRow, iCol) = saRecovery;
                    
                end
            end
            
            % returning baseline
            baseline.mapMean = mapBaseline;
            baseline.dataMean = nanmean(mapBaseline(:));
            
            % returning amplitude
            amplitude.mapMean = mapAmplitude;
            amplitude.dataMean = nanmean(mapAmplitude(:));
            
            % returning duration
            duration.mapMean = mapDuration;
            duration.dataMean = nanmean(mapDuration(:));
            
            % returning data on activation maps
            activationMaps.mapMean = mapActivationTimes;
            
            % returning data on recovery maps
            recoveryMaps.mapMean = mapRecoveryTimes;
            
        end
        
        %%%%% POSTPROCESSING %%%%%
        function alternans = getAlternans(data, varargin)
            % Extracts alternans quantity for a given feature (e.g. amplitude or duration).
            % It can either work on a vector of numbers, giving alternans between odd/even values,
            % or on a stack of multiple wave passes (giving spatial
            % alternans map over odd/even slices).
            %
            % IN:
            % data - either a vector of numbers or a stack of spatial maps
            % of the feature on which alternans is to be measured (e.g.
            % amplitude.maps produced by COSMAS.analyseRegularPacing).
            %
            % varargin - an optional parameter which may determine the method for alternans estimation. In
            % both, average for even and odd values/slices is computed. Then,
            % 'largerToSmaller' (default) measures ratio of larger to smaller, and
            % 'sMAPE' does abs(odd-even)/(odd+even).
            %
            % OUT:
            % alternans - if data is a vector, this gives a single number,
            % if data is a stack of spatial maps, it produces a single
            % spatial map.
            method = 'largerToSmaller';
            if (~isempty(varargin))
                method = varargin{1};
            end
            assert(ismember(method, {'largerToSmaller', 'sMAPE'}), 'The parameter specifying the method must be either largerToSmaller or sMAPE.');
            
            if (length(size(data))==2) % if we have a vector of numbers, we convert it to a stack so we can process everything the same way
                data = COSMAS.traceToStack(data);
            end
            
            meanOdd = nanmean(data(:,:,1:2:end), 3);
            meanEven = nanmean(data(:,:,2:2:end), 3);
            
            maxMap = max(meanOdd,meanEven);
            minMap = min(meanOdd,meanEven);
            if (isequal(method, 'largerToSmaller'))
                alternans = maxMap ./ minMap;
            elseif (isequal(method, 'sMAPE'))
                alternans = abs(meanOdd - meanEven)./(meanOdd + meanEven);
            end
        end
        
        function cv = getCV(activationMap, XY, varargin)
            % Analyse conduction velocity (CV) between pair (or pairs) of points.
            % By default, this returns the CV in pixels per frame, but can also
            % provide it in cm/s.
            %
            % IN:
            % activationMap - a single activation map.
            %
            % XY - a matrix of size n-by-4 encoding pairs of points between which CV is measured using the provided activation map.
            % Each row corresponds to an origin and target point (columns
            % are: rowFrom, columnFrom, rowTo, columnTo).
            %
            % varargin may optionally contain a two-numbers parameter allowing specification of spatial (how many mm is a single pixel side) and temporal resolution (in frames per second). If not given, the output is in pixels/frame, otherwise in cm/s.
            %
            % OUT:
            % cv -  a vector of conduction velocities, one per row of XY.
            nRows = size(XY, 1);
            cv = zeros(nRows, 1);
            
            for iRow = 1:nRows
                fromRow = XY(iRow, 1);
                fromCol = XY(iRow, 2);
                toRow = XY(iRow, 3);
                toCol = XY(iRow, 4);
                distDiff = sqrt((fromRow - toRow)^2 + (fromCol - toCol)^2);
                timeDiff = activationMap(toRow, toCol) - activationMap(fromRow, fromCol);
            end
            
            cv = distDiff/timeDiff;
            
            % Potential rescaling to cm/s
            if (~isempty(varargin))
                assert(length(varargin{1})==2, 'the scaling vector must be a 2-by-1 or 1-by-2 vetor specifying how many mm is a single pixel and what is the temporal resolution (fps)');
                distancePerPixel = varargin{1}(1);
                fps = varargin{1}(2);
                cv = (cv * distancePerPixel*fps)/10; % /10 is to convert mm/ms to cm/s
            end
        end
        
        function xyuv = getLocalCV(activationMap, baylyNeighbourhood, varargin)
            % Performs local estimation of CV using Bayly's method (doi: 10.1109/10.668746), producing a vector field and optionally its plot.
            %
            % IN: 
            % activationMap - activation map from which a vector field is
            % obtained.
            %
            % baylyNeighbourhood - distance around a point that is
            % considered when fitting the Bayly polynomial.
            %
            % varargin{1} - if defined, gives maximum length of a velocity
            % vector (the longer ones are discarded).
            %
            % varargin{2} - if defined, contains the index of figure in
            % which the CV field is drawn. If not defined, no figure is
            % produced.
            %
            % varargin{3} - if defined, the output path of storage of
            % varargin{2}.
            %
            % OUT:
            % xyuv - a n-by-4 matrix, where n is number of pixels and columns correspond
            % to x,y,u,v:  x,y gives indices of row and column, with u,v
            % corresponding to dx,dy. Mind that this is in row/column
            % coordinates - if plotting via quiver (in standard x-y
            % coordinates), this needs to be shuffled slightly, see the code
            % at the end of the function.
            maxDistance = Inf;
            if (length(varargin)>=1)
                maxDistance = varargin{1}; % maximum arrow length that is not discarded
            end
            
            figureNumber = [];
            if (length(varargin)>=2)
                figureNumber = varargin{2}; % maximum arrow length that is not discarded
            end
            
            foutName = [];
            if (length(varargin)>=3)
                foutName = varargin{3}; % maximum arrow length that is not discarded
            end
            
            %% Bayly CV estimation
            nRows = size(activationMap, 1);
            nCols = size(activationMap, 2);
            [m,n] = ndgrid(1:nRows, 1:nCols);
            Z = [m(:), n(:)];
            activationTimes = activationMap(:); % same order of activation times as in Z.
            
            xyuv = zeros(size(Z,1), 4);
            
            for iPoint = 1:size(Z,1)
                % for each point, find points nearby.
                thisPoint = Z(iPoint,:);
                distances = pdist2(thisPoint, Z);
                whereNeighbours = find((distances>=0) & (distances<=baylyNeighbourhood));
                locationsNeighbours = Z(whereNeighbours, :);
                neighbourActivationTimes = activationTimes(whereNeighbours);
                try
                    sf = fit([locationsNeighbours], neighbourActivationTimes, 'poly22');
                catch
                    xyuv(iPoint, :) = [thisPoint(1), thisPoint(2), NaN, NaN]; % Massive vector that will be rejected
                    continue;
                end

                coeffs = coeffvalues(sf); %const, x, y, x2, xy, y2 ;
                % IN PAPER: f, d, e, a, c, b
                x = thisPoint(1);
                y = thisPoint(2);
                dx = coeffs(2)+2*coeffs(4)*x + coeffs(5)*y;  % x
                dy = coeffs(3)+2*coeffs(6)*y + coeffs(5)*x;  % y
                xyuv(iPoint, :) = [x, y, dx/(dx*dx+dy*dy), dy/(dx*dx+dy*dy)];
            end
            
            %% We get rid of dx,dy which are too long
            arrowLengths = sqrt(xyuv(:,3).*xyuv(:,3) + xyuv(:,4).*xyuv(:,4));
            xyuv(arrowLengths > maxDistance, 3) = NaN;
            xyuv(arrowLengths > maxDistance, 4) = NaN;
            
            %% Optional plotting of the quiver. Given that axes are
            % different between x/y (and for Matlab row/columns, [1,1] is
            % top left, while for common set of axes it is bottom left, so
            % we do complement for the 2nd parameter in quiver and we flip
            % the sign of the fourth parameter in quiver)
            
            if (~isempty(figureNumber))
                figure(figureNumber); 
                quiver(xyuv(:,2), 16-xyuv(:,1), xyuv(:,4), -xyuv(:,3), 0, 'k');
                
                if (~isempty(foutName))
                    saveas(gcf, foutName);
                end
            end
        end
        
        function plotActivationMap(activationMap, varargin)
            % This function plots a contour map of activation.
            %
            % IN:
            % activationMap - a single activation map from which contours are
            % obtained.
            %
            % varagin{1} - optionally, the index of figure used for this purpose may be given - if not
            % specified, a new figure is opened and used.
            %
            % varargin{2} - optional parameter defining the
            % 'levels' parameter of @contourf (basically, a contour line is
            % drawn each varargin{2} ms).
            figureNumber = [];
            if (~isempty(varargin))
                figureNumber = varargin{1};
                assert(isnumeric(figureNumber), 'the second parameter (figure number) must be a positive integer');
                assert(figureNumber > 0, 'the second parameter (figure number) must be a positive integer');
            end
            
            levelGranularity = 1;
            if (length(varargin)>=2)
                levelGranularity = varargin{2};
               
            end
            
            if isempty(figureNumber)
                figure;
            else
                figure(figureNumber);
            end
            contourf(flipud(activationMap - min(min(activationMap))), 'levels', levelGranularity, 'ShowText', 'on');
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
        end
        
        %%%%% HELPER FUNCTIONS %%%%%
        function imStackOut = applyMask(imStackIn, mask)
            % This function may be used to apply a binary mask to each frame of a 3D stack, setting pixels-to-be-discarded to NaN.
            % This may be useful to get rid of empty space around the image
            % of the heart, space around a Petri dish for cultures, etc.
            %
            % IN:
            % imStackIn - a 3D stack representing a recording.
            %
            % mask - a binary mask of the same size as a single frame of
            % imStackIn. Where mask==0, the corresponding pixels are set to
            % NaN.
            %
            % OUT:
            % imStackOut - imStackIn where zero elements in mask are set to
            % NaN for each frame.
            
            imStackOut = imStackIn;
            
            for iFrame = 1:size(imStackIn, 3)
                frame = imStackIn(:,:,iFrame);
                frame(mask==0) = NaN;
                imStackOut(:,:,iFrame) = frame;
            end
        end
        
        function binnedImage = binningVec(img, logfactor)
            % A function which carries out spatial binning for an image
            % (replacing each bin with the average of the pixels found in
            % it).
            %
            % IN:
            % img - an image to be spatially binned.
            %
            % logfactor - base-two logarithm of the binning factor - this must be a
            % positive integer (i.e., logfactor of 1 leads to 2-by-2
            % binning, 2 to 4-by4, 3 to 8-by-8, etc.). Note that while this
            % helper function requires a logarithm of the binning factor,
            % the data-reading functions readFolder and readTifStack take
            % the binning factor directly, computing its logarithm
            % internally, so that the user doesn't have to take care of
            % that.
            %
            % OUT:
            % binnedImage - the source image after spatial binning.
            binnedImage = img;
            for j = 1:logfactor % How many binnings
                oddRows = double(binnedImage(1:2:end,:)); % A submatrix of blkimg of odd columns
                evenRows = double(binnedImage(2:2:end,:)); % A submatrix of blkimg of even columns
                rowBinned = (oddRows + evenRows) / 2;
                oddCols = rowBinned(:, 1:2:end); % And similar thing done for rows
                evenCols = rowBinned(:, 2:2:end);
                binnedImage = (oddCols + evenCols)/2; % img binned by 2^j
            end
        end
        
        function vectorMinima = combGetMinima(signalTrace, bcl, varargin)
            % A function searching for minima in traces from cardiac preparations with known activation pattern.
            % The function can be naturally also used to extract signal maxima when
            % the source trace is inverted.
            %
            % IN:
            % signalTrace - a vector containing signal, such as calcium transients or action potentials
            % (e.g., intensity of a pixel in optical mapping, or an electrophysiological recording).
            %
            % bcl - the basic cycle length (number of frames between two activations).
            %
            % varargin{1} - The first optional parameter is the refinementWidth parameter for comb algorithm
            % (in ms/samples - the radius of local search around comb
            % teeth). Default is 10.
            %
            % varargin{2} - The second optional parameter is the custom comb that is used instead of the
            % regularly placed one. It should be an increasing vector of
            % numbers, where the first element determines the last possible
            % position of the first minimum to be searched for, and the
            % subsequent elements give further indices of minima. E.g.
            % using [30 180, 330, 400] means that the algorithm will search
            % for 3 minima that are 150 frames apart, and one that is
            % further 100 frames after the last previous one, and the first
            % minimum is to be placed between frame 1 and frame 30. See
            % sampleScript_5_S1S2 for an example of custom comb. If this parameter is defined, then
            % the parameter ‘bcl’ is not used within the code (it still has to be provided,
            % but it can be any number).
            %
            %
            % OUT:
            % vectorMinima - the vector of local minima in the signal, which are
            % approximately bcl ms (or samples) apart.
            
            %% Default parameter initialization and processing of extra inputs.
            refinementWidth = 10;
            customComb = [];
            
            if (numel(varargin) == 1) % 1st extra argument is the width of local search
                if (~isempty(varargin{1}))
                    refinementWidth = varargin{1};
                end
            elseif (numel(varargin) == 2)
                if (~isempty(varargin{1}))
                    refinementWidth = varargin{1};
                end
                if (~isempty(varargin{2}))
                    customComb = varargin{2};
                end
                
            end
            
            
            %% Comb positioning
            nFrames = length(signalTrace);
            
            if (isempty(customComb)) % if no custom comb is given, we use a regular one based on bcl.
                candidateMinimaMean = inf*ones(bcl-1, 1);
                for iStart = 1:bcl
                    candidateMinima = iStart:bcl:nFrames;
                    candidateMinimaMean(iStart) = mean(signalTrace(candidateMinima));
                end
            else % when custom comb is provided, we slide it up to the first element
                candidateMinimaMean = inf * ones(customComb(1), 1);
                combNoOffset = customComb - customComb(1); % comb with first element subtracted so it starts at 0.
                for iStart = 1:customComb(1)
                    candidateMinima = iStart + combNoOffset;
                    candidateMinimaMean(iStart) = mean(signalTrace(candidateMinima));
                end
            end
            
            % The start leading to minimal average value is used as the
            % indicator of minima.
            [~, bestMinima] = min(candidateMinimaMean);
            if (isempty(customComb))
                vectorMinima =  bestMinima:bcl:nFrames;
            else
                vectorMinima = bestMinima + combNoOffset;
            end
            
            %% Comb refinement
            % The vector of teeth (~near-minima) is refined by looking refinementWidth to the
            % left and right, taking the actual minimum there.
            vectorMinimaRefined = zeros(size(vectorMinima));
            
            for iMinimum = 1:length(vectorMinima)
                leftBound = max(vectorMinima(iMinimum) - refinementWidth, 1);
                rightBound = min(vectorMinima(iMinimum) + refinementWidth, length(signalTrace));
                
                [~, actualMinimum] = min(signalTrace(leftBound:rightBound));
                vectorMinimaRefined(iMinimum) = leftBound + actualMinimum - 1;
            end
            vectorMinima = vectorMinimaRefined;
            
        end
        
        function [apd, timeRecovery] = getAPD(time, activationSignal, level, varargin)
            % The function returns the duration of an action potential (or calcium transient)
            % at the desired level of repolarization.
            % The function assumes the action potential/calcium transient to be upward-pointing (i.e. peak activation
            % is more positive than resting value) - if
            % used on downward-pointing signal, this needs to be inverted first.
            %
            % IN:
            % time - the time vector for the activation trace.
            %
            % activationSignal - the trace of a single action potential or
            % calcium transient (or any similar signal).
            %
            % level - the level of recovery at which the duration is to be
            % measured. This is scaled between 0 and 1 (i.e., for APD80,
            % the value 0.8 is to be used).
            %
            % varargin{1} - a string encoding the method of baseline
            % estimation; see the documentation of COSMAS.analyseRegularPacing,
            % the 'parameters.baselineDefinition' parameter.
            %
            % varargin{2} - a string encoding the method of object
            % detection/selection; see the documentation of
            % COSMAS.analyseRegularPacing, the 'parameters.objectDetection'
            % parameter.
            %
            % varargin{3} - when 'augmented' objectDetection is used, this
            % parameter is used to pass the time of peak upstroke time
            % within the single activation provided in 'activationSignal'.
            % When more objects are above the threshold determined by
            % 'level', the one closest to this parameter is used.
            %
            % OUT:
            %
            % apd - the duration of the action potential/calcium transient
            % at the given level of recovery. The function uses interpolation to get sub-frame resolution.
            %
            % timeRecovery - the time of the end of the action
            % potential/calcium transient at the given repolarization
            % level. The function uses interpolation to get sub-frame resolution.
            
            baselineDefinition = 'first';
            objectDetection = 'first';
            if (length(varargin)>=1)
                baselineDefinition = varargin{1};
            end
            if length(varargin)>=2
                objectDetection = varargin{2};
            end
            if length(varargin)==3
                peakDiffActivation = varargin{3};
            end
            
            
            baseline = activationSignal(1);
            if (isequal(baselineDefinition, 'firstlast'))
                baseline = mean(activationSignal([1,end]));
            end
            thresh = (baseline + (1-level)*(max(activationSignal)-baseline));
            s = bwconncomp(activationSignal > thresh);
            
            % Now, the longest object is found among those above threshold.
            % It's almost always the first one, but if there is a lot of
            % noise (or pepper noise), it's not necessarily the case.
            
            %% we switch between taking the largest and first found object
            if (s.NumObjects > 1) && isequal(objectDetection, 'largest')
                lengths = cellfun(@length, s.PixelIdxList);
                [~, whereMax] = max(lengths);
                interval = s.PixelIdxList{whereMax};
            elseif (s.NumObjects > 1) && (strcmp(objectDetection, 'augmented'))
                assert(~isempty(peakDiffActivation), 'when augmented object detection used, peakDiffActivation parameter (3rd varargin) must be specified');
                distances = zeros(s.NumObjects, 1);
                for iObject = 1:s.NumObjects
                    pix = s.PixelIdxList{iObject};
                    distances(iObject) = min(abs(pix-peakDiffActivation));
                end
                [~, whereMin] = min(distances);
                interval = s.PixelIdxList{whereMin};
            else
                interval = s.PixelIdxList{1};
            end
            
            % linear interpolation of times when the signal crosses
            % threshold at the activation and recovery is carried out
            tStartInterval = time(interval(1));
            tEndInterval = time(interval(end));
            try % when activation duration reaches until the end of recording, it is impossible to find the subsequent element in time vector - this pretty much always indicates a problem with the data, e.g. when starting value is 1, then peak goes to 100, and then repolarization stabilizes at 50 (can happen with movement artefacts or very noisy data) - then it artefactually takes almost the whole trace as action potential. Same for the activation starting so high that the first frame is already considered "active" - then the signal is potentially incomplete. We don't want such artefactual cases to return numerical values, so NaN is returned.
                tBeforeStartInterval = time(interval(1)-1);
                tAfterEndInterval = time(interval(end)+1);
                
            catch
                apd = NaN;
                timeRecovery = NaN;
                return;
            end
            % processing interval before activation start
            vStartInterval = activationSignal(interval(1));
            vBeforeStartInterval = activationSignal(interval(1) - 1);
            howFarThreshBeforeStart = (vStartInterval-thresh) / (vStartInterval-vBeforeStartInterval);
            interpAddingStart = (tStartInterval - tBeforeStartInterval) * howFarThreshBeforeStart;
            % processing interval after recovery
            vEndInterval = activationSignal(interval(end));
            vAfterEndInterval = activationSignal(max(interval) + 1);
            howFarThreshAfterEnd = (thresh-vEndInterval) / (vAfterEndInterval-vEndInterval);
            interpAddingEnd = (tAfterEndInterval - tEndInterval) * howFarThreshAfterEnd;
            apd = time(interval(end)) - time(interval(1)) + interpAddingStart + interpAddingEnd ;
            
            timeRecovery = time(interval(end)) + interpAddingEnd;
        end
        
        function halfActivationTime = getHalfActivationTime(activationSignal, varargin)
            % The function returns the time at "half-activation" of an action potential
            % or calcium transient; taken with regards to amplitude (i.e. it corresponds
            % to the time when APD50 would be measured from when applied to action potentials).
            % It uses linear interpolation to get this half-activation time
            % "accurately" (as in not just the previous/subsequent frame).
            % The function assumes the action potential/calcium transient to be upward-pointing
            % (i.e. peak activation is more positive than resting value) - if used
            % on downward-pointing signal, this needs to be inverted first.
            %
            % IN:
            % activationSignal - the trace of a single action potential or
            % calcium transient (or any similar signal).
            %
            % varargin{1} - a string encoding the method of baseline
            % estimation; see the documentation of COSMAS.analyseRegularPacing,
            % the 'parameters.baselineDefinition' parameter.
            %
            % varargin{2} - a string encoding the method of object
            % detection/selection; see the documentation of
            % COSMAS.analyseRegularPacing, the 'parameters.objectDetection'
            % parameter.
            %
            % varargin{3} - when 'augmented' objectDetection is used, this
            % parameter is used to pass the time of peak upstroke time
            % within the single activation provided in 'activationSignal'.
            % When more objects are above the threshold determined by
            % 'level', the one closest to this parameter is used.
            %
            % OUT:
            % halfActivationTime - The time of half-activation of the signal
            % (i.e. time when 50% of the signal amplitude is achieved).
            % The function uses interpolation to get sub-frame resolution.
            
            baselineDefinition = 'first';
            objectDetection = 'first';
            peakDiffActivation = [];
            if (~isempty(varargin))
                baselineDefinition = varargin{1};
            end
            if length(varargin)>=2
                objectDetection = varargin{2};
            end
            if length(varargin)==3
                peakDiffActivation = varargin{3};
            end
            
            
            try
                baseline = activationSignal(1);
                if (isequal(baselineDefinition, 'firstlast'))
                    baseline = mean(activationSignal([1,end]));
                end
                peak = min(activationSignal); % peak potential. This is minimum given ap signal heads down
                threshold = baseline - 0.5 * (baseline - peak);
                s = bwconncomp(activationSignal < threshold);
                
                if (s.NumObjects > 1) && strcmp(objectDetection, 'largest') % picking largest object as the activation to be measured
                    lengths = cellfun(@length, s.PixelIdxList);
                    [~, whereMax] = max(lengths);
                    apFrames = s.PixelIdxList{whereMax};
                elseif (s.NumObjects > 1) && (strcmp(objectDetection, 'augmented'))
                    assert(~isempty(peakDiffActivation), 'when augmented object detection used, peakDiffActivation parameter (3rd varargin) must be specified');
                    distances = zeros(s.NumObjects, 1);
                    for iObject = 1:s.NumObjects
                        pix = s.PixelIdxList{iObject};
                        distances(iObject) = min(abs(pix-peakDiffActivation));
                    end
                    [~, whereMin] = min(distances);
                    apFrames = s.PixelIdxList{whereMin};
                else
                    apFrames = s.PixelIdxList{1};
                end
                
                afterHalfActivationTime = apFrames(1);
                preHalfActivationTime = apFrames(1) - 1;
                traceAfter = activationSignal(afterHalfActivationTime);
                tracePre = activationSignal(preHalfActivationTime);
                halfActivationTime = preHalfActivationTime + ((threshold - tracePre)/(traceAfter - tracePre)) * (afterHalfActivationTime - preHalfActivationTime);% the last parentheses are 1 ms at the moment - written like this just for future generality.
            catch
                disp('issue in getting half activation')
                halfActivationTime = NaN;
            end
        end
        
        function imOut = getImageForMask(imStack)
            % This function gives a high-contrast average image of the stack (averaged over
            % time). This can be used to draw a mask (manually or automatically).
            %
            % IN:
            % imStack - A stack representing the recording. The input can be of type
            % uint8, uint16, or double (that one is expected to be scaled between
            % 0 and 1 - if not, please convert the stack to the uint types).
            %
            % OUT:
            % imOut - A maximum-contrast image of the time-average of the
            % stack.
            avgIm = mean(imStack, 3);
            avgIm = avgIm - min(avgIm(:)); % subtracting baseline
            maxIm = max(avgIm(:));
            if (isequal(class(avgIm), 'double'))
                imOut = avgIm * 1/maxIm;
            elseif (isequal(class(avgIm), 'uint8'))
                imOut = avgIm * 255/maxIm;
            elseif (isequal(class(avgIm), 'uint16'))
                imOut = avgIm * 65535/maxIm;
            else
                error('Stack type must be double, uin8, uint16.');
            end
        end
        
        function [saBaseline, saAmplitude, saDuration, saActivation, saRecovery] = processSingleActivation(activationTrace, spikesPointDown, baselineDefinition,durationLevel, objectDetection,recordingClock, timeStart, peakDiffActivation, verbose)
            % Processes a single activation (action potential or calcium
            % transient), returning its properties/features - this is
            % mainly a helper function for COSMAS.analyseRegularPacing or
            % COSMAS.analyseSinglePass, and the parameters it uses are the
            % same as there.
            
            % Sorting where data will be stored within the
            % coordinates given by recordingClock.
            if (spikesPointDown)
                [~, wherePeak] = min(activationTrace);
            else
                [~, wherePeak] = max(activationTrace);
            end
            
            % Try-catch, because when signal is poor,
            % getHalfActivation may fail, and in that case,
            % we don't really want the other values either
            try
                % Baseline - first index
                bas = activationTrace(1);
                if (isequal(baselineDefinition, 'firstlast'))
                    bas = mean(activationTrace([1,end]));
                end
                saBaseline = bas;
                
                % Amplitude
                saAmplitude = abs(activationTrace(wherePeak) - bas); % abs takes care of spikes pointing up/down
                % Duration
                if (spikesPointDown)
                    [saDuration, saRecovery] = COSMAS.getAPD(1:length(activationTrace), -activationTrace, durationLevel, baselineDefinition, objectDetection, peakDiffActivation);
                else
                    [saDuration, saRecovery] = COSMAS.getAPD(1:length(activationTrace), activationTrace, durationLevel, baselineDefinition, objectDetection, peakDiffActivation);
                end
                
                % Activation - first extracted in local coordinates, converted to global
                % coordinates, and then it is taken relatively
                % to the clock
                if (spikesPointDown)
                    activationLocal = COSMAS.getHalfActivationTime(activationTrace, baselineDefinition, objectDetection, peakDiffActivation);
                else
                    activationLocal = COSMAS.getHalfActivationTime(-activationTrace, baselineDefinition, objectDetection, peakDiffActivation);
                end
                activationGlobal = timeStart + activationLocal - 1;
                whereLastSmaller = find(recordingClock < activationGlobal, 1, 'last');
                
                if (isempty(whereLastSmaller) || (whereLastSmaller >= length(recordingClock)))
                    error('There is a problem with data, possibly with the rate of activation, or there may be an arrhythmia');
                end
                saActivation = activationGlobal - recordingClock(whereLastSmaller); % computing activation relative to preceding global clock
                
                % handling recovery similarly to activation, matching it to
                % the global clock. If APD extraction failed (because the
                % cell started or ended the trace activated), NaN is
                % returned.
                %
                if (~isnan(saRecovery))
                    recoveryGlobal = timeStart + saRecovery - 1;
                    saRecovery = recoveryGlobal - recordingClock(whereLastSmaller); % here we reuse the clock information from when getting relative activation time
                else
                    saRecovery = NaN;
                end
            catch me
                if (verbose)
                    disp('problem processing a trace of a pixel');
                    disp(me.message)
                end
                saBaseline = NaN;
                saAmplitude = NaN;
                saDuration = NaN;
                saActivation = NaN;
                saRecovery = NaN;
            end
        end
        
        function [traceSmoothed, traceSmoothedDetrended] = smoothTrace(tracePixel, filterWidth, varargin)
            % This function smooths (~denoises) a trace using Savitzky-Golay filter of
            % a given width. It may also remove drift/trend in data using
            % polynomial fitting.
            %
            % IN:
            % tracePixel - the trace to be smoothed.
            %
            % filterWidth - the width of the Savitzky-Golay filter (should
            % be an odd number - if an even one is provided, one is added
            % to it).
            %
            % varargin{1} - the optional parameter may contain the order of the polynomial
            % used to detrend the signal (useful for removing drift etc.).
            % Unlike standard methods which subtract the polynomial making
            % the signal approximately zero-centered on the y-axis, here we
            % re-add the average of the original trace to the zero-centered
            % trace, so that the information on signal baseline is not
            % lost.
            %
            % OUT:
            % traceSmoothed - tracePixel after smoothing.
            %
            % traceSmoothedDetrended - tracePixel after smoothing and
            % baseline subtraction.
            
            tracePixel = tracePixel(:); % making sure this is a column vector
            %% detrending
            toDetrend = ~isempty(varargin);
            
            if toDetrend
                detrendOrder = varargin{1};
                assert(isnumeric(detrendOrder), 'Order of baseline subtraction must be numeric');
                assert(detrendOrder>0, 'Order of baseline subtraction must be larger than 0');
                
                % below, we fit a polynomial to the signal, which acts as a
                % baseline. Mu shenanigans serve to improve the numerical
                % stability (see documentation of polyfit for mu).
                [p,s,mu] = polyfit([1:length(tracePixel)]', tracePixel, detrendOrder); % mean subtraction is just so that the problem isn't ill-conditioned
                polyBaseline = polyval(p, [1:length(tracePixel)]', [], mu);
                traceDeTrend = tracePixel - polyBaseline + mean(tracePixel); %we re-add mean of tracePixel so that we don't lose information on total signal baseline
            end
            
            if (mod(filterWidth, 2) == 0)
                filterWidth = filterWidth + 1;
            end
            
            if (filterWidth > 4) % if filter width is <= 4, it's not possible to carry out filtering because the filter width must be more than filter order. In this way, we can also forbid filtering (by setting smoothingParameter to 1)
                traceSmoothed = sgolayfilt(tracePixel,4,filterWidth);
                traceSmoothedDetrended = traceSmoothed; % just so that this is defined as an output parameter if no detrending is used.
                if toDetrend
                    traceSmoothedDetrended = sgolayfilt(traceDeTrend,4,filterWidth);
                end
            else
                traceSmoothed = tracePixel;
                traceSmoothedDetrended = traceSmoothed;
                if toDetrend
                    traceSmoothedDetrended = traceDeTrend;
                end
            end
            
            
        end
        
        function imStack = traceToStack(signalTrace)
            % A minor helper function which converts a vector (either n-by-1 or 1-by-n)
            % to a “fake” 3D stack of 1-by-1-by-n – in this way, it is the same format
            % as a stack representing a video (albeit with 1x1 resolution) and it can
            % be fed to analyseRegularPacing. One can therefore use the functionality
            % of COSMAS with regards to baseline, amplitude, or duration features even
            % on single traces (coming e.g. from electrophysiological recordings).
            %
            % IN:
            % signalTrace - the trace to be converted into a "stack".
            %
            % OUT:
            % imStack - the resulting stack.
            imStack = zeros(1, 1, length(signalTrace));
            imStack(:) = signalTrace;
        end
        %%%%% NATSORT ETC. %%%%%
        function [X,ndx,dbg] = natsortfiles(X,rgx,varargin)
            % Alphanumeric / Natural-Order sort of a cell array of filename/filepath strings (1xN char).
            %
            % (c) 2014-2019 Stephen Cobeldick
            %
            % Alphanumeric sort of a cell array of filenames or filepaths: sorts by
            % character order and also by the values of any numbers that are within
            % the names. Filenames, file-extensions, and directories (if supplied)
            % are split apart and are sorted separately: this ensures that shorter
            % filenames sort before longer ones (i.e. thus giving a dictionary sort).
            %
            %%% Example:
            % P = 'C:\SomeDir\SubDir';
            % S = dir(fullfile(P,'*.txt'));
            % C = natsortfiles({S.name});
            % for k = 1:numel(C)
            %     fullfile(P,C{k})
            % end
            %
            %%% Syntax:
            %  Y = natsortfiles(X)
            %  Y = natsortfiles(X,rgx)
            %  Y = natsortfiles(X,rgx,<options>)
            % [Y,ndx,dbg] = natsortfiles(X,...)
            %
            % To sort all of the strings in a cell array use NATSORT (File Exchange 34464).
            % To sort the rows of a cell array of strings use NATSORTROWS (File Exchange 47433).
            %
            %% File Dependency %%
            %
            % NATSORTFILES requires the function NATSORT (File Exchange 34464). The optional
            % arguments <options> are passed directly to NATSORT. See NATSORT for case
            % sensitivity, sort direction, numeric substring matching, and other options.
            %
            %% Explanation %%
            %
            % Using SORT on filenames will sort any of char(0:45), including the printing
            % characters ' !"#$%&''()*+,-', before the file extension separator character '.'.
            % Therefore this function splits the name and extension and sorts them separately.
            %
            % Similarly the file separator character within filepaths can cause longer
            % directory names to sort before shorter ones, as char(0:46)<'/' and
            % char(0:91)<'\'. Check this example to see why this matters:
            %
            % >> X = {'A1\B', 'A+/B', 'A\B', 'A=/B', 'A/B'};
            % >> sort(X)
            % ans =   'A+/B'  'A/B'  'A1\B'  'A=/B'  'A\B'
            % >> natsortfiles(X)
            % ans =   'A\B'  'A/B'  'A1\B'  'A+/B'  'A=/B'
            %
            % NATSORTFILES splits filepaths at each file separator character and sorts
            % every level of the directory hierarchy separately, ensuring that shorter
            % directory names sort before longer, regardless of the characters in the names.
            %
            %% Examples %%
            %
            % >> A = {'a2.txt', 'a10.txt', 'a1.txt'};
            % >> sort(A)
            % ans = 'a1.txt'  'a10.txt'  'a2.txt'
            % >> natsortfiles(A)
            % ans = 'a1.txt'  'a2.txt'  'a10.txt'
            %
            % >> B = {'test_new.m'; 'test-old.m'; 'test.m'};
            % >> sort(B) % Note '-' sorts before '.':
            % ans =
            %    'test-old.m'
            %    'test.m'
            %    'test_new.m'
            % >> natsortfiles(B) % Shorter names before longer (dictionary sort):
            % ans =
            %    'test.m'
            %    'test-old.m'
            %    'test_new.m'
            %
            % >> C = {'test2.m'; 'test10-old.m'; 'test.m'; 'test10.m'; 'test1.m'};
            % >> sort(C) % Wrong numeric order:
            % ans =
            %    'test.m'
            %    'test1.m'
            %    'test10-old.m'
            %    'test10.m'
            %    'test2.m'
            % >> natsortfiles(C) % Shorter names before longer:
            % ans =
            %    'test.m'
            %    'test1.m'
            %    'test2.m'
            %    'test10.m'
            %    'test10-old.m'
            %
            %%% Directory Names:
            % >> D = {'A2-old\test.m';'A10\test.m';'A2\test.m';'A1archive.zip';'A1\test.m'};
            % >> sort(D) % Wrong numeric order, and '-' sorts before '\':
            % ans =
            %    'A10\test.m'
            %    'A1\test.m'
            %    'A1archive.zip'
            %    'A2-old\test.m'
            %    'A2\test.m'
            % >> natsortfiles(D) % Shorter names before longer (dictionary sort):
            % ans =
            %    'A1archive.zip'
            %    'A1\test.m'
            %    'A2\test.m'
            %    'A2-old\test.m'
            %    'A10\test.m'
            %
            %% Input and Output Arguments %%
            %
            %%% Inputs (*=default):
            % X   = CellArrayOfCharRowVectors, with filenames or filepaths to be sorted.
            % rgx = Regular expression to match number substrings, '\d+'*
            %     = [] uses the default regular expression, which matches integers.
            % <options> can be supplied in any order and are passed directly to NATSORT.
            %
            %%% Outputs:
            % Y   = CellArrayOfCharRowVectors, filenames of <X> sorted into natural-order.
            % ndx = NumericMatrix, same size as <X>. Indices such that Y = X(ndx).
            % dbg = CellVectorOfCellArrays, size 1xMAX(2+NumberOfDirectoryLevels).
            %       Each cell contains the debug cell array for directory names, filenames,
            %       and file extensions. Helps debug the regular expression. See NATSORT.
            %
            % See also SORT NATSORT NATSORTROWS DIR FILEPARTS FULLFILE NEXTNAME CELLSTR REGEXP IREGEXP SSCANF
            
            %% Input Wrangling %%
            %
            assert(iscell(X),'First input <X> must be a cell array.')
            tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
            assert(all(tmp(:)),'First input <X> must be a cell array of strings (1xN character).')
            %
            if nargin>1
                varargin = [{rgx},varargin];
            end
            %
            %% Split and Sort File Names/Paths %%
            %
            % Split full filepaths into file [path,name,extension]:
            [pth,fnm,ext] = cellfun(@fileparts,X(:),'UniformOutput',false);
            % Split path into {dir,subdir,subsubdir,...}:
            pth = regexp(pth,'[^/\\]+','match'); % either / or \ as filesep.
            len = cellfun('length',pth);
            num = max(len);
            vec = cell(numel(len),1);
            %
            % Natural-order sort of the file extensions and filenames:
            if isempty(num)
                ndx = [];
                ids = [];
                dbg = {};
            elseif nargout<3 % faster:
                [~,ndx] = COSMAS.natsort(ext,varargin{:});
                [~,ids] = COSMAS.natsort(fnm(ndx),varargin{:});
            else % for debugging:
                [~,ndx,dbg{num+2}] = COSMAS.natsort(ext,varargin{:});
                [~,ids,tmp] = COSMAS.natsort(fnm(ndx),varargin{:});
                [~,idd] = sort(ndx);
                dbg{num+1} = tmp(idd,:);
            end
            ndx = ndx(ids);
            %
            % Natural-order sort of the directory names:
            for k = num:-1:1
                idx = len>=k;
                vec(:) = {''};
                vec(idx) = cellfun(@(c)c(k),pth(idx));
                if nargout<3 % faster:
                    [~,ids] = COSMAS.natsort(vec(ndx),varargin{:});
                else % for debugging:
                    [~,ids,tmp] = COSMAS.natsort(vec(ndx),varargin{:});
                    [~,idd] = sort(ndx);
                    dbg{k} = tmp(idd,:);
                end
                ndx = ndx(ids);
            end
            %
            % Return the sorted array and indices:
            ndx = reshape(ndx,size(X));
            X = X(ndx);
            %
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortfiles
        
        function [X,ndx,dbg] = natsort(X,rgx,varargin)
            % Alphanumeric / Natural-Order sort the strings in a cell array of strings (1xN char).
            %
            % (c) 2012-2019 Stephen Cobeldick
            %
            % Alphanumeric sort a cell array of strings: sorts by character order and
            % also by the values of any number substrings. Default: match all integer
            % number substrings and perform a case-insensitive ascending sort.
            %
            %%% Example:
            % >> X = {'x2', 'x10', 'x1'};
            % >> sort(X)
            % ans =   'x1'  'x10'  'x2'
            % >> natsort(X)
            % ans =   'x1'  'x2'  'x10'
            %
            %%% Syntax:
            %  Y = natsort(X)
            %  Y = natsort(X,rgx)
            %  Y = natsort(X,rgx,<options>)
            % [Y,ndx,dbg] = natsort(X,...)
            %
            % To sort filenames or filepaths use NATSORTFILES (FEX 47434).
            % To sort the rows of a cell array of strings use NATSORTROWS (FEX 47433).
            %
            %% Number Substrings %%
            %
            % By default consecutive digit characters are interpreted as an integer.
            % Specifying the optional regular expression pattern allows the numbers to
            % include a +/- sign, decimal digits, exponent E-notation, quantifiers,
            % or look-around matching. For information on defining regular expressions:
            % http://www.mathworks.com/help/matlab/matlab_prog/regular-expressions.html
            %
            % The number substrings are parsed by SSCANF into numeric values, using
            % either the *default format '%f' or the user-supplied format specifier.
            %
            % This table shows examples of regular expression patterns for some common
            % notations and ways of writing numbers, with suitable SSCANF formats:
            %
            % Regular       | Number Substring | Number Substring              | SSCANF
            % Expression:   | Match Examples:  | Match Description:            | Format Specifier:
            % ==============|==================|===============================|==================
            % *         \d+ | 0, 123, 4, 56789 | unsigned integer              | %f  %i  %u  %lu
            % --------------|------------------|-------------------------------|------------------
            %      [-+]?\d+ | +1, 23, -45, 678 | integer with optional +/- sign| %f  %i  %d  %ld
            % --------------|------------------|-------------------------------|------------------
            %     \d+\.?\d* | 012, 3.45, 678.9 | integer or decimal            | %f
            % (\d+|Inf|NaN) | 123, 4, NaN, Inf | integer, Inf, or NaN          | %f
            %  \d+\.\d+e\d+ | 0.123e4, 5.67e08 | exponential notation          | %f
            % --------------|------------------|-------------------------------|------------------
            %  0[0-7]+      | 012, 03456, 0700 | octal notation & prefix       | %o  %i
            %   [0-7]+      |  12,  3456,  700 | octal notation                | %o
            % --------------|------------------|-------------------------------|------------------
            %  0X[0-9A-F]+  | 0X0, 0X3E7, 0XFF | hexadecimal notation & prefix | %x  %i
            %    [0-9A-F]+  |   0,   3E7,   FF | hexadecimal notation          | %x
            % --------------|------------------|-------------------------------|------------------
            %  0B[01]+      | 0B1, 0B101, 0B10 | binary notation & prefix      | %b   (not SSCANF)
            %    [01]+      |   1,   101,   10 | binary notation               | %b   (not SSCANF)
            % --------------|------------------|-------------------------------|------------------
            %
            %% Debugging Output Array %%
            %
            % The third output is a cell array <dbg>, to check if the numbers have
            % been matched by the regular expression <rgx> and converted to numeric
            % by the SSCANF format. The rows of <dbg> are linearly indexed from <X>,
            % even columns contain numbers, odd columns contain split substrings:
            %
            % >> [~,~,dbg] = natsort(X)
            % dbg =
            %    'x'    [ 2]
            %    'x'    [10]
            %    'x'    [ 1]
            %
            %% Examples %%
            %
            %%% Multiple integers (e.g. release version numbers):
            % >> A = {'v10.6', 'v9.10', 'v9.5', 'v10.10', 'v9.10.20', 'v9.10.8'};
            % >> sort(A)
            % ans =   'v10.10'  'v10.6'  'v9.10'  'v9.10.20'  'v9.10.8'  'v9.5'
            % >> natsort(A)
            % ans =   'v9.5'  'v9.10'  'v9.10.8'  'v9.10.20'  'v10.6'  'v10.10'
            %
            %%% Integer, decimal, NaN, or Inf numbers, possibly with +/- signs:
            % >> B = {'test+NaN', 'test11.5', 'test-1.4', 'test', 'test-Inf', 'test+0.3'};
            % >> sort(B)
            % ans =   'test' 'test+0.3' 'test+NaN' 'test-1.4' 'test-Inf' 'test11.5'
            % >> natsort(B, '[-+]?(NaN|Inf|\d+\.?\d*)')
            % ans =   'test' 'test-Inf' 'test-1.4' 'test+0.3' 'test11.5' 'test+NaN'
            %
            %%% Integer or decimal numbers, possibly with an exponent:
            % >> C = {'0.56e007', '', '43E-2', '10000', '9.8'};
            % >> sort(C)
            % ans =   ''  '0.56e007'  '10000'  '43E-2'  '9.8'
            % >> natsort(C, '\d+\.?\d*([eE][-+]?\d+)?')
            % ans =   ''  '43E-2'  '9.8'  '10000'  '0.56e007'
            %
            %%% Hexadecimal numbers (with '0X' prefix):
            % >> D = {'a0X7C4z', 'a0X5z', 'a0X18z', 'a0XFz'};
            % >> sort(D)
            % ans =   'a0X18z'  'a0X5z'  'a0X7C4z'  'a0XFz'
            % >> natsort(D, '0X[0-9A-F]+', '%i')
            % ans =   'a0X5z'  'a0XFz'  'a0X18z'  'a0X7C4z'
            %
            %%% Binary numbers:
            % >> E = {'a11111000100z', 'a101z', 'a000000000011000z', 'a1111z'};
            % >> sort(E)
            % ans =   'a000000000011000z'  'a101z'  'a11111000100z'  'a1111z'
            % >> natsort(E, '[01]+', '%b')
            % ans =   'a101z'  'a1111z'  'a000000000011000z'  'a11111000100z'
            %
            %%% Case sensitivity:
            % >> F = {'a2', 'A20', 'A1', 'a10', 'A2', 'a1'};
            % >> natsort(F, [], 'ignorecase') % default
            % ans =   'A1'  'a1'  'a2'  'A2'  'a10'  'A20'
            % >> natsort(F, [], 'matchcase')
            % ans =   'A1'  'A2'  'A20'  'a1'  'a2'  'a10'
            %
            %%% Sort order:
            % >> G = {'2', 'a', '', '3', 'B', '1'};
            % >> natsort(G, [], 'ascend') % default
            % ans =   ''   '1'  '2'  '3'  'a'  'B'
            % >> natsort(G, [], 'descend')
            % ans =   'B'  'a'  '3'  '2'  '1'  ''
            % >> natsort(G, [], 'num<char') % default
            % ans =   ''   '1'  '2'  '3'  'a'  'B'
            % >> natsort(G, [], 'char<num')
            % ans =   ''   'a'  'B'  '1'  '2'  '3'
            %
            %%% UINT64 numbers (with full precision):
            % >> natsort({'a18446744073709551615z', 'a18446744073709551614z'}, [], '%lu')
            % ans =       'a18446744073709551614z'  'a18446744073709551615z'
            %
            %% Input and Output Arguments %%
            %
            %%% Inputs (*==default):
            % X   = CellArrayOfCharRowVectors, to be sorted into natural-order.
            % rgx = Regular expression to match number substrings, '\d+'*
            %     = [] uses the default regular expression, which matches integers.
            % <options> can be entered in any order, as many as required:
            %     = Sort direction: 'descend'/'ascend'*
            %     = NaN/number order: 'NaN<num'/'num<NaN'*
            %     = Character/number order: 'char<num'/'num<char'*
            %     = Character case handling: 'matchcase'/'ignorecase'*
            %     = SSCANF number conversion format, e.g.: '%f'*, '%x', '%li', '%b', etc.
            %
            %%% Outputs:
            % Y   = CellArrayOfCharRowVectors, <X> sorted into natural-order.
            % ndx = NumericArray, such that Y = X(ndx). The same size as <X>.
            % dbg = CellArray of the parsed characters and number values.
            %       Each row is one input char vector, linear-indexed from <X>.
            %
            % See also SORT NATSORTFILES NATSORTROWS CELLSTR REGEXP IREGEXP SSCANF
            
            %% Input Wrangling %%
            %
            assert(iscell(X),'First input <X> must be a cell array.')
            tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
            assert(all(tmp(:)),'First input <X> must be a cell array of char row vectors (1xN char).')
            %
            if nargin<2 || isnumeric(rgx)&&isempty(rgx)
                rgx = '\d+';
            else
                assert(ischar(rgx)&&ndims(rgx)<3&&size(rgx,1)==1,...
                    'Second input <rgx> must be a regular expression (char row vector).') %#ok<ISMAT>
            end
            %
            % Optional arguments:
            tmp = cellfun('isclass',varargin,'char') & cellfun('size',varargin,1)<2 & cellfun('ndims',varargin)<3;
            assert(all(tmp(:)),'All optional arguments must be char row vectors (1xN char).')
            % Character case:
            ccm = strcmpi(varargin,'matchcase');
            ccx = strcmpi(varargin,'ignorecase')|ccm;
            % Sort direction:
            sdd = strcmpi(varargin,'descend');
            sdx = strcmpi(varargin,'ascend')|sdd;
            % Char/num order:
            chb = strcmpi(varargin,'char<num');
            chx = strcmpi(varargin,'num<char')|chb;
            % NaN/num order:
            nab = strcmpi(varargin,'NaN<num');
            nax = strcmpi(varargin,'num<NaN')|nab;
            % SSCANF format:
            sfx = ~cellfun('isempty',regexp(varargin,'^%([bdiuoxfeg]|l[diuox])$'));
            %
            COSMAS.nsAssert(1,varargin,sdx,'Sort direction')
            COSMAS.nsAssert(1,varargin,chx,'Char<->num')
            COSMAS.nsAssert(1,varargin,nax,'NaN<->num')
            COSMAS.nsAssert(1,varargin,sfx,'SSCANF format')
            COSMAS.nsAssert(0,varargin,~(ccx|sdx|chx|nax|sfx))
            %
            % SSCANF format:
            if nnz(sfx)
                fmt = varargin{sfx};
                if strcmpi(fmt,'%b')
                    cls = 'double';
                else
                    cls = class(sscanf('0',fmt));
                end
            else
                fmt = '%f';
                cls = 'double';
            end
            %
            %% Identify Numbers %%
            %
            [mat,spl] = regexpi(X(:),rgx,'match','split',varargin{ccx});
            %
            % Determine lengths:
            nmx = numel(X);
            nmn = cellfun('length',mat);
            nms = cellfun('length',spl);
            mxs = max(nms);
            %
            % Preallocate arrays:
            bon = bsxfun(@le,1:mxs,nmn).';
            bos = bsxfun(@le,1:mxs,nms).';
            arn = zeros(mxs,nmx,cls);
            ars =  cell(mxs,nmx);
            ars(:) = {''};
            ars(bos) = [spl{:}];
            %
            %% Convert Numbers to Numeric %%
            %
            if nmx
                tmp = [mat{:}];
                if strcmp(fmt,'%b')
                    tmp = regexprep(tmp,'^0[Bb]','');
                    vec = cellfun(@(s)sum(pow2(s-'0',numel(s)-1:-1:0)),tmp);
                else
                    vec = sscanf(sprintf(' %s',tmp{:}),fmt);
                end
                assert(numel(vec)==numel(tmp),'The %s format must return one value for each input number.',fmt)
            else
                vec = [];
            end
            %
            %% Debugging Array %%
            %
            if nmx && nargout>2
                dbg = cell(mxs,nmx);
                dbg(:) = {''};
                dbg(bon) = num2cell(vec);
                dbg = reshape(permute(cat(3,ars,dbg),[3,1,2]),[],nmx).';
                idf = [find(~all(cellfun('isempty',dbg),1),1,'last'),1];
                dbg = dbg(:,1:idf(1));
            else
                dbg = {};
            end
            %
            %% Sort Columns %%
            %
            if ~any(ccm) % ignorecase
                ars = lower(ars);
            end
            %
            if nmx && any(chb) % char<num
                boe = ~cellfun('isempty',ars(bon));
                for k = reshape(find(bon),1,[])
                    ars{k}(end+1) = char(65535);
                end
                [idr,idc] = find(bon);
                idn = sub2ind(size(bon),boe(:)+idr(:),idc(:));
                bon(:) = false;
                bon(idn) = true;
                arn(idn) = vec;
                bon(isnan(arn)) = ~any(nab);
                ndx = 1:nmx;
                if any(sdd) % descending
                    for k = mxs:-1:1
                        [~,idx] = sort(COSMAS.nsGroup(ars(k,ndx)),'descend');
                        ndx = ndx(idx);
                        [~,idx] = sort(arn(k,ndx),'descend');
                        ndx = ndx(idx);
                        [~,idx] = sort(bon(k,ndx),'descend');
                        ndx = ndx(idx);
                    end
                else % ascending
                    for k = mxs:-1:1
                        [~,idx] = sort(ars(k,ndx));
                        ndx = ndx(idx);
                        [~,idx] = sort(arn(k,ndx),'ascend');
                        ndx = ndx(idx);
                        [~,idx] = sort(bon(k,ndx),'ascend');
                        ndx = ndx(idx);
                    end
                end
            else % num<char
                arn(bon) = vec;
                bon(isnan(arn)) = ~any(nab);
                if any(sdd) % descending
                    [~,ndx] = sort(COSMAS.nsGroup(ars(mxs,:)),'descend');
                    for k = mxs-1:-1:1
                        [~,idx] = sort(arn(k,ndx),'descend');
                        ndx = ndx(idx);
                        [~,idx] = sort(bon(k,ndx),'descend');
                        ndx = ndx(idx);
                        [~,idx] = sort(COSMAS.nsGroup(ars(k,ndx)),'descend');
                        ndx = ndx(idx);
                    end
                else % ascending
                    [~,ndx] = sort(ars(mxs,:));
                    for k = mxs-1:-1:1
                        [~,idx] = sort(arn(k,ndx),'ascend');
                        ndx = ndx(idx);
                        [~,idx] = sort(bon(k,ndx),'ascend');
                        ndx = ndx(idx);
                        [~,idx] = sort(ars(k,ndx));
                        ndx = ndx(idx);
                    end
                end
            end
            %
            ndx  = reshape(ndx,size(X));
            X = X(ndx);
            %
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsort
        function nsAssert(val,inp,idx,varargin)
            % Throw an error if an option is overspecified.
            if nnz(idx)>val
                tmp = {'Unknown input arguments',' option may only be specified once. Provided inputs'};
                error('%s:%s',[varargin{:},tmp{1+val}],sprintf('\n''%s''',inp{idx}))
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsAssert
        function grp = nsGroup(vec)
            % Groups of a cell array of strings, equivalent to [~,~,grp]=unique(vec);
            [vec,idx] = sort(vec);
            grp = cumsum([true,~strcmp(vec(1:end-1),vec(2:end))]);
            grp(idx) = grp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsGroup
    end
end

