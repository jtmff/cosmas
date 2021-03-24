import os
import numpy as np
from scipy import signal
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from imageio import imread
from natsort import natsorted
from skimage import measure
import math
import tifffile
import code

class COSMAS:
##### READING DATA #####
    def readFolder(self,folderName,varargin=None):
        # a helper function returning a stack of images from folder with
        # image sequence. Uses natsort to get natural ordering (this is useful when there are
        # files 1.png, 2.png,...,10.png,11.png, which would get sorted as 1,10,11,2, in purely alphabetic sorting, which is not the right order of frames)
        #
        # IN:
        # folderName is the folder containing a sequence of images
        # corresponding to frames.
        #
        # varargin can contain a structure of parameters (called e.g. readerParameters), with the
        # following fields:
        # 1) readerParameters.extension - extension of images in the
        # given folder.
        # 2) readerParameters.fromTo - a 2-by-1 vector [from to],
        # determining from which to which frame is the recording to be
        # processed.
        # 3) readerParameters.binningFactor - spatial binning parameter
        # (must be a power of two) - binningFactor x binningFactor
        # pixels are aggregated into a single one.
        #
        # OUT:
        # imStack - the resulting 3D stack, where rows and columns
        # correspond to frame rows and columns, and 3rd dimension to
        # index of image in the sequence (i.e. time).
        #
        # avgTrace - a trace corresponding to spatially averaged stack.


        ## Setting defaults, reading parameters, getting filenames to be read.
        extension = '.tif'
        fromTo = []
        binningFactor = 1
        if varargin:
            readerParameters = varargin
            if hasattr(readerParameters,'extension'):
                extension = readerParameters.extension
            if hasattr(readerParameters,'fromTo'):
                fromTo = readerParameters.fromTo
            if hasattr(readerParameters,'binningFactor'):
                binningFactor = readerParameters.binningFactor
        fnames = []
        for file in os.listdir(folderName):
            if file.endswith(extension):
                fnames.append(file)
        fnames = natsorted(fnames)

        if not fromTo:
            fromTo = [0, len(fnames)]

        ## Checking parameter correctness
        assert fromTo[0]>=1 & fromTo[1] <= len(fnames),'Boundaries from-to the frame are incorrect. Fewer images are available than the boundaries ask for.'
        assert round(np.log2(binningFactor)) == np.log2(binningFactor), 'The binning factor must be a power of 2.'

        ## Reading stack
        im = imread(folderName + '/' + fnames[0])
        nFrames = fromTo[1] - fromTo[0]

        imStack = np.zeros((int(np.shape(im)[0]/binningFactor),int(np.shape(im)[1]/binningFactor), nFrames))
        avgTrace = np.zeros((nFrames,1))

        for iFrame in range(fromTo[0], fromTo[1]):
            img = imread(folderName+'/'+fnames[iFrame])
            if (binningFactor != 1):
                img = self.binningVec(img,np.log2(binningFactor))
            imStack[:,:,iFrame-fromTo[0]] = img
            imStack = imStack.astype('uint16')
            avgTrace[iFrame-fromTo[0]] = np.nanmean(img[:])
        varargout = avgTrace
        return imStack, varargout

    def readTifStack(self, fname, varargin=None):
        # returning a stack of images from a single tif stack.
        # fname is the path to the file (including filename and extension).
        #
        # IN:
        # varargin can contain a structure of parameters (called readerParameters), with the
        # following fields:
        # 1) readerParameters.fromTo - a 2-by-1 vector [from to],
        # determining from which to which frame is the recording to be
        # processed
        # 2) readerParameters.binningFactor - spatial binning parameter
        # (must be a power of two) - binningFactor x binningFactor
        # pixels are aggregated into a single one.
        #
        # OUT:
        # imStack - the resulting 3D stack.
        #
        # avgTrace - a trace corresponding to spatially averaged stack.



        imgs = tifffile.TiffFile(fname)
        nFramesTotal = len(imgs.pages)

        fromTo = []
        binningFactor = 1
        if varargin is not None:
            readerParameters = varargin
        if hasattr(readerParameters,'fromTo'):
            fromTo = readerParameters.fromTo
        if hasattr(readerParameters,'binningFactor'):
            binningFactor = readerParameters.binningFactor
        if not fromTo:
            fromTo = [1, nFramesTotal]

        ## Checking parameter correctness
        assert fromTo[0]>=1 & fromTo[1] <= nFramesTotal, 'Boundaries from-to the frame are incorrect. Fewer images are available than the boundaries ask for.'
        assert round(np.log2(binningFactor)) == np.log2(binningFactor), 'The binning factor must be a power of 2.'

        ## Reading stack
        nFrames = fromTo[1] - fromTo[0]
        imStack = np.zeros((int(imgs.pages[0].shape[0]/binningFactor), int(imgs.pages[0].shape[1]/binningFactor),nFrames))
        avgTrace = np.zeros((nFrames,1))

        iFrom = fromTo[0]
        iTo = fromTo[1]
        imgs = imgs.asarray()
        for iFrame in range(iFrom, iTo):
            img = imgs[iFrame]
            if (binningFactor != 1):
                img = self.binningVec(img, np.log2(binningFactor))
            imStack[:,:,iFrame - iFrom] = img
            #imStack = imStack.astype('uint16')
            avgTrace[iFrame-iFrom] = self.nanmean(img[:])
        varargout = avgTrace
        return imStack, varargout


##### ANALYSING DATA #####
    def analyseRegularPacing(self, imageStack, bcl, parameters):
        # A function for processing stacks corresponding to
        # recordings with multiple passes of a wave.
        #
        # IN
        # imageStack - a stack representing the recording.
        #
        # bcl - basic cycle length of the recording (in frames).
        #
        # parameters - a structure with parameters:
        #
        # parameters.baselineDefinition - whether baseline of an
        # activation (e.g. calcium transient) is taken as the first
        # element in each segmented activation ('first'), or as the
        # average of the first and last element ('firstlast').
        #
        # parameters.baselineSubtractionOrder - order of polynomial
        # subtraction of signal baseline. Use -1 if no baseline
        # subtraction is to be done (or just leave the parameter field
        # undefined).
        #
        # parameters.durationLevel - level at which duration is
        # extracted, scaled to 0-1. E.g., use 0.8 for APD80.
        #
        # parameters.objectDetection - when multiple activations are discovered in
        # a single signal segment (while only one can be true calcium
        # transient/action potential), this parameter determines how
        # the correct activation is detected (hopefully :)). 'first' -
        # first object found in the segment. 'largest' - the largest
        # object is picked (with most frames). 'augmented' -
        # information on derivative of the signal is used, see the
        # publication. For voltage mapping with multiple wave passes,
        # we recommend using 'augmented' (which is not great for calcium,
        # given that there are no sharp upstrokes). Otherwise, 'largest'
        # is usually more robust than 'first'.
        #
        # parameters.smoothingParameter - width of Savitzky-Golay
        # filtering for signal smoothing. It should be an odd number
        # (if an even number is given, 1 will be added). Given that 4th
        # order smoothing is used, this parameter has to be at least 5
        # when provided (when smaller, no filtering is done).
        #
        # parameters.spikesPointDown - if true, signal activation
        # manifests as reduction in signal intensity (e.g. some voltage
        # dyes), i.e., action potentials "point down". default =
        # false.
        #
        # parameters.verbose - if true, the code reports when there is
        # a problem with segmentation and/or processing of a pixel
        # trace (which can happen when a pixel contains only noise, for
        # example).
        #
        # parameters.waveProcessing – if 'perbeat', each wave is processed separately,
        # and the output structures contain a map for each complete wave pass.
        # If 'hybrid', a recording clock is used to chop the recording into sub-stacks, w
        # hich are then averaged, and a single-wave processing is applied to this subsequently.
        # The value 'hybrid' is good for very noisy recordings and activation mapping,
        # but is not suggested to be used for APD mapping or amplitude measurements.
        # If ‘hybrid’ is used, parameter.objectDetection should be set to ‘largest’.
        # The default is ‘perbeat’.
        #
        # OUT:
        # baseline - a structure describing signal baseline (e.g. bases
        # of calcium transients
        #
        # amplitude - a structure describing signal amplitude
        #
        # duration - a structure describing signal duration (e.g. APD)
        #
        # activationMaps - a structure describing activation pattern in
        # the recording (relative to recording clock; if you want to
        # have minimum activation in 0, just subtract minimum of these maps)
        #
        # recoveryMaps - a structure describing recovery pattern (e.g.
        # the time when APD80 is reached). Relative to the same clock
        # as activationMaps
        #
        # recordingClock - the recording clock determining global
        # synchronization of single segmented activations
        #
        # the structures baseline, amplitude, and duration have the
        # following fields (shown for duration):
        #    duration.maps - a 3D stack where each slice corresponds to
        #    a map of the feature in a single wave pass
        #
        #    duration.mapMean - the average of the previous maps
        #
        #    duration.data - a vector of spatial averages of maps in duration.maps
        #
        #    duration.dataMean = nanmean(duration.data) - the mean of
        #    the field data (i.e., this gives one number summarizing
        #    the whole recording)
        #
        #    duration.mapMeanEven - mean map for even beats (useful for
        #    inspection of alternans).
        #
        #    duration.mapMeanOdd - mean map for odd beats
        #
        #    duration.dataMeanEven - spatial averages of even maps of
        #    the feature
        #
        #    duration.dataMeanOdd - spatial averages of odd maps of the
        #    feature
        #
        # the structure activationMaps and recoveryMaps have the same
        # fields, except the ones starting with 'data' (there is not
        # much point in spatially averaged activation).

        ## Reading and verifying parameters
        spikesPointDown = False
        if hasattr(parameters,'spikesPointDown'):
            spikesPointDown = parameters.spikesPointDown

        baselineDefinition = 'first'
        if hasattr(parameters, 'baselineDefinition'):
            baselineDefinition = parameters.baselineDefinition
        def ismember(A,B):
            return [np.sum(a == B) for a in A ]
        assert ismember(baselineDefinition, ['first','firstlast']), 'baselineDefinition must be either first or firstlast'

        objectDetection = 'largest'
        if hasattr(parameters,'objectDetection'):
            objectDetection = parameters.objectDetection
        assert ismember(objectDetection,['augmented','first','largest']), 'objectDetection must be either first or largest'

        baselineSubtractionOrder = -1
        if hasattr(parameters,'baselineSubtractionOrder'):
            baselineSubtractionOrder = parameters.baselineSubtractionOrder

        smoothingParameter = 11
        if hasattr(parameters, 'smoothingParameter'):
            smoothingParameter = parameters.smoothingParameter

        durationLevel = 0.75
        if hasattr(parameters, 'durationLevel'):
            durationLevel = parameters.durationLevel

        verbose = False
        if hasattr(parameters, 'verbose'):
            verbose = parameters.verbose
        assert isinstance(verbose,bool), 'verbose parameter must be either true or false'

        waveProcessing = 'perbeat'
        if hasattr(parameters, 'waveProcessing'):
            waveProcessing = parameters.waveProcessing
        assert ismember(waveProcessing,['perbeat', 'hybrid']), 'waveProcessing must be perbeat or hybrid'

        customComb = []
        if hasattr(parameters, 'customComb'):
            customComb = parameters.customComb

        nRows = np.shape(imageStack)[0]
        nCols = np.shape(imageStack)[1]

        # diastolicDetector = 'comb'
        # if hasattr(parameters,'diastolicDetector'):
        #     diastolicDetector = parameters.diastolicDetector
        # assert ismember(diastolicDetector,['comb','prePeak']), 'diastolicDetector must be either first or firstlast'

        ## Getting global signal, from which we extract the recording 'timer' that separates action potentials/calcium transients
        # Unlike usual, we extract locations of peaks p1,p2,..., and then build a
        # timer that starts at (p1+p2)/2, so peaks happen roughly at
        # midpoint between the points in the timer
        avgTrace = np.squeeze(self.nanmean(self.nanmean(imageStack,0),0))
        if (baselineSubtractionOrder > 0):
            avgTrace = self.smoothTrace(avgTrace, smoothingParameter, baselineSubtractionOrder)[1]
        else:
            avgTrace = self.smoothTrace(avgTrace, smoothingParameter)

        if (spikesPointDown):
            maxActivations = self.combGetMinima(avgTrace, bcl, [[], customComb])
        else:
            maxActivations = self.combGetMinima(-avgTrace, bcl, [[], customComb])

        recordingClock = np.zeros(len(maxActivations)-1)
        for i in range(0,len(maxActivations)-1):
            recordingClock[i] = self.nanmean([maxActivations[i], maxActivations[i+1]])

        if customComb:
            recordingClock = list(recordingClock)
            recordingClock.append(np.shape(imageStack)[2]-1)
            recordingClock.insert(0, 0)
            recordingClock = np.array(recordingClock)
        if (waveProcessing == 'hybrid'): # hybrid processing, where the timer is used to chop the recording to smaller parts and it is then averaged and processed using singleWavePass
            recordingClock = np.arange(recordingClock[0],len(avgTrace),bcl)

            # Now, we process all traces in the stack with smoothing
            # and baseline subtraction, then chopping it to bcl-sized
            # chunks and averaging them, before processing it as a
            # single wave
            for iRow in range(0, nRows):
                for iCol in range(0, nCols):
                    # extracting and smoothing the trace
                    pixelTrace = np.squeeze(imageStack[iRow,iCol,:])
                    if (sum(np.isnan(pixelTrace))>0):
                        continue

                    if (baselineSubtractionOrder > 0):
                        traceSmoothed = self.smoothTrace(pixelTrace, smoothingParameter, baselineSubtractionOrder)[1]
                    else:
                        traceSmoothed = self.smoothTrace(pixelTrace, smoothingParameter)
                    imageStack[iRow, iCol, :] = traceSmoothed
            stackAvg = np.zeros((nRows, nCols, bcl))
            for iPass in range(0, len(recordingClock)-1):
                stackAvg = stackAvg + imageStack[:,:,np.arange(recordingClock[iPass], recordingClock[iPass+1], 1).astype(int)]
                # traces[iPass] = avgTrace[np.arange(recordingClock[iPass],recordingClock[iPass+1], 1).astype(int)]
            stackAvg = stackAvg / (len(recordingClock) - 1)
            baseline, amplitude, duration, activationMaps, recoveryMaps = self.analyseSinglePass(stackAvg, parameters)
            return baseline, amplitude, duration, activationMaps, recoveryMaps, recordingClock

        else: # otherwise we do standard processing, beat per beat
            mapsBaseline = np.full((nRows, nCols, len(recordingClock)-1),np.nan)
            mapsAmplitude = np.full((nRows, nCols, len(recordingClock)-1),np.nan)
            mapsDuration = np.full((nRows, nCols, len(recordingClock)-1),np.nan)
            mapsActivationTimes = np.full((nRows, nCols, len(recordingClock)-1),np.nan)
            mapsRecoveryTimes = np.full((nRows, nCols, len(recordingClock)-1),np.nan)

            # for each trace, we segment it using comb, and use local
            # maxima to assign the found values using recordingClock. We
            # don't want to use recordingClock or anything like that itself
            # for segmentation of spikes/transients, as that does not give
            # fine-enough information (e.g. in discordant alternans)
            for iRow in range(0, nRows):
                for iCol in range(0, nCols):
                    pixelTrace = np.squeeze(imageStack[iRow,iCol,:])
                    if (sum(np.isnan(pixelTrace))>0):
                        continue

                    if (baselineSubtractionOrder > 0):
                        traceSmoothed = self.smoothTrace(pixelTrace, smoothingParameter, baselineSubtractionOrder)[1]
                    else:
                        traceSmoothed = self.smoothTrace(pixelTrace, smoothingParameter)

                    # For the trace, we extract properties between all its diastoles
                    if (spikesPointDown):
                        diastoles = self.combGetMinima(-1*traceSmoothed, bcl, [[], customComb])
                    else:
                        diastoles = self.combGetMinima(traceSmoothed, bcl, [[], customComb])
                    # we also find minima/maxima of dv/dt that serve
                    # 'augmented' object detection (finding objects nearest peak)
                    # diff (signal).
                    diffSignal = self.smoothTrace(np.diff(traceSmoothed), smoothingParameter)
                    if (spikesPointDown):
                        peakDiffs = self.combGetMinima(diffSignal, bcl, [[], customComb])
                    else:
                        peakDiffs = self.combGetMinima(-1*diffSignal, bcl, [[], customComb])

                    for iStart in range(0, len(diastoles)-1):
                        timeStart = int(diastoles[iStart])
                        timeEnd = int(diastoles[iStart + 1])

                        iCenter = (timeStart + timeEnd)/2.0 # converting the peak within single activation transient to the global temporal coordinates.
                        iBin = sum(iCenter > recordingClock)# after how many elements of the recording clock does the location come?
                        if (iBin<1) or (iBin >= len(recordingClock)):
                            continue
                        iBin = iBin - 1

                        peakDiffs = np.array(peakDiffs).astype(int)
                        peakDiffActivation = peakDiffs[np.where((peakDiffs>=timeStart) & (peakDiffs <= timeEnd))] - timeStart + 1 # Take time of peak activation relative to time start

                        # Single activation is extracted and processed
                        activationTrace = traceSmoothed[timeStart:timeEnd+1]

                        saBaseline, saAmplitude, saDuration, saActivation, saRecovery = self.processSingleActivation(activationTrace, spikesPointDown, baselineDefinition, durationLevel, objectDetection, recordingClock, timeStart, peakDiffActivation, verbose)


                        mapsBaseline[iRow,iCol,iBin] = saBaseline
                        mapsAmplitude[iRow,iCol,iBin] = saAmplitude
                        mapsDuration[iRow,iCol,iBin] = saDuration
                        mapsActivationTimes[iRow,iCol,iBin] = saActivation
                        mapsRecoveryTimes[iRow,iCol,iBin] = saRecovery

            ## Now we remove slices with empty entries in activation - this
            # refers to some pixels having not enough information (e.g. not
            # enough time at the end to contain a full action potential/CaT), so
            # these slices are discarded.

            nZerosInSlice = sum(sum(mapsActivationTimes==0))
            mapsBaseline = mapsBaseline[:,:,nZerosInSlice == 0]
            mapsAmplitude = mapsAmplitude[:, :, nZerosInSlice == 0]
            mapsDuration = mapsDuration[:,:,nZerosInSlice == 0]
            mapsActivationTimes = mapsActivationTimes[:,:,nZerosInSlice == 0]

            class output:
                maps = 0.0
                mapMean = 0.0
                data = 0.0
                dataMean = 0.0
                mapMeanEven = 0.0
                mapMeanOdd = 0.0
                dataMeanEven = 0.0
                dataMeanOdd = 0.0

            baseline = output()
            amplitude  = output()
            duration = output()
            activationMaps = output()
            recoveryMaps = output()

            def saveOutputInStructure(x, outV):
                x.maps = outV
                x.mapMean = self.nanmean(outV,2)
                x.data = np.array([np.squeeze(self.nanmean(self.nanmean(outV,0),0))]).transpose()
                x.dataMean = [self.nanmean(x.data,0)]
                if (np.shape(outV)[2]>1):
                    end = len(outV[0,0,:])
                    x.mapMeanEven = self.nanmean(outV[:,:,1:end:2],2)
                    x.mapMeanOdd = self.nanmean(outV[:,:,0:end:2],2)
                    end = len(x.data)
                    x.dataMeanEven = np.squeeze(self.nanmean(x.data[1:end:2]))
                    x.dataMeanOdd = np.squeeze(self.nanmean(x.data[0:end:2]))

            # Returning baseline
            saveOutputInStructure(baseline, mapsBaseline)
            saveOutputInStructure(amplitude, mapsAmplitude)
            saveOutputInStructure(duration, mapsDuration)
            saveOutputInStructure(activationMaps, mapsActivationTimes)
            saveOutputInStructure(recoveryMaps, mapsRecoveryTimes)

            # returning clock separating activations. This is filtered to
            # remove elements corresponding to incompletely filled maps
            # (i.e. usually this would be just the last map, when it's not
            # completely filled because some pixels don't have a
            # long-enough signal there)
            clockFiltered = recordingClock
            sliceDetected = np.where(nZerosInSlice > 0)[0]
            if (sliceDetected):
                clockFiltered[sliceDetected + 1] = []
            varargout = clockFiltered
            return baseline, amplitude, duration, activationMaps, recoveryMaps,varargout

    def analyseSinglePass(self, imageStack, parameters):
        # A function for processing stacks corresponding to
        # recordings with single pass of a wave. Please see the
        # documentation of SCUMS.analyseRegularPacing for the
        # description of inputs/outputs.
        #
        # In the output structures, the only fields are mapMean and
        # dataMean. This may look illogical (there is just one map per
        # the recording, so why would one carry out averaging over slices,
        # which is essentially identity?), but
        # it's based on our experience that SCUMS users mainly use the
        # outputs of this function in a similar way as they'd use
        # mapMean/dataMean returned by analyseRegularPacing -
        # therefore, calling the outputs the same means there is less
        # code rewriting needed when one switches between
        # analyseRegularPacing and analyseSinglePass.

        spikesPointDown = False # if spikes point down instead of up. This is relevant for finding diastole, and for extracting baseline of the signal.
        if hasattr(parameters,'spikesPointDown'):
            spikesPointDown = parameters.spikesPointDown

        baselineDefinition = 'first'
        if hasattr(parameters, 'baselineDefinition'):
            baselineDefinition = parameters.baselineDefinition

        def ismember(A,B):
            return [np.sum(a == B) for a in A ]
        assert ismember(baselineDefinition, ['first','firstlast']), 'baselineDefinition must be either first or firstlast'

        objectDetection = 'largest'
        if hasattr(parameters,'objectDetection'):
            objectDetection = parameters.objectDetection
        assert ismember(objectDetection,['augmented','first','largest']), 'objectDetection must be either first or largest'

        baselineSubtractionOrder = -1
        if hasattr(parameters,'baselineSubtractionOrder'):
            baselineSubtractionOrder = parameters.baselineSubtractionOrder

        smoothingParameter = 11
        if hasattr(parameters, 'smoothingParameter'):
            smoothingParameter = parameters.smoothingParameter

        durationLevel = 0.75
        if hasattr(parameters, 'durationLevel'):
            durationLevel = parameters.durationLevel

        verbose = False
        if hasattr(parameters, 'verbose'):
            verbose = parameters.verbose
        assert isinstance(verbose,bool), 'verbose parameter must be either true or false'

        nRows = np.shape(imageStack)[0]
        nCols = np.shape(imageStack)[1]
        mapsBaseline = np.full((nRows, nCols),np.nan)
        mapsAmplitude = np.full((nRows, nCols),np.nan)
        mapsDuration = np.full((nRows, nCols),np.nan)
        mapsActivationTimes = np.full((nRows, nCols),np.nan)
        mapsRecoveryTimes = np.full((nRows, nCols),np.nan)
        for iRow in range(0, nRows):
            for iCol in range(0, nCols):
                pixelTrace = np.squeeze(imageStack[iRow,iCol,:])
                # we smooth the trace, but do not perform baseline
                # subtraction - that can cause a huge mess in
                # single-wave traces
                traceSmoothed = self.smoothTrace(pixelTrace, smoothingParameter)
                if (sum(np.isnan(pixelTrace))>0):
                    continue
                timeStart = 1
                diffSignal = np.diff(traceSmoothed) + 1
                if (spikesPointDown):
                    peakDiffActivation = np.where(diffSignal==diffSignal.min())[0][-1]
                else:
                    peakDiffActivation = np.where(-diffSignal==-diffSignal.min())[0][-1]
                recordingClock = np.array([0, np.inf]) # this is a fairly dummy value, making sure that the results are considered to belong to the first wave (out of 1).
                # single activation is extracted and processed
                saBaseline, saAmplitude, saDuration, saActivation, saRecovery = self.processSingleActivation(traceSmoothed, spikesPointDown, baselineDefinition, durationLevel, objectDetection, recordingClock, timeStart, peakDiffActivation, verbose)
                mapsBaseline[iRow,iCol] = saBaseline
                mapsAmplitude[iRow,iCol] = saAmplitude
                mapsDuration[iRow,iCol] = saDuration
                mapsActivationTimes[iRow,iCol] = saActivation
                mapsRecoveryTimes[iRow,iCol] = saRecovery
        class output:
            mapMean = 0.0
            dataMean = 0.0
        baseline = output()
        amplitude = output()
        duration = output()
        activationMaps = output()
        recoveryMaps = output()
        baseline.mapMean = mapsBaseline
        baseline.dataMean = [self.nanmean(mapsBaseline[:],0)]
        amplitude.mapMean = mapsAmplitude
        amplitude.dataMean = [self.nanmean(mapsAmplitude[:],0)]
        duration.mapMean = mapsDuration
        duration.dataMean = [self.nanmean(mapsDuration[:],0)]
        activationMaps.mapMean = mapsActivationTimes
        recoveryMaps.mapMean = mapsRecoveryTimes
        return baseline, amplitude, duration, activationMaps, recoveryMaps

##### POSTPROCESSING #####
    def getAlternans(self, data, varargin):
        # Extracts alternans quantity for a given feature (e.g. amplitude or duration).
        # It can either work on a vector of numbers, giving alternans between odd/even values,
        # or on a stack of multiple wave passes (giving spatial
        # alternans map over odd/even slices).
        #
        # IN:
        # data - either a vector of numbers or a stack of spatial maps
        # of the feature on which alternans is to be measured (e.g.
        # amplitude.maps produced by SCUMS.analyseRegularPacing).
        #
        # varargin - an optional parameter which may determine the method for alternans estimation. In
        # both, average for even and odd values/slices is computed. Then,
        # 'largerToSmaller' (default) measures ratio of larger to smaller, and
        # 'sMAPE' does abs(odd-even)/(odd+even).
        #
        # OUT:
        # alternans - if data is a vector, this gives a single number,
        # if data is a stack of spatial maps, it produces a single
        # spatial map.
        method = 'largerToSmaller'
        if varargin:
            method = varargin
        assert method in ['largerToSmaller','sMAPE'], 'The parameter specifying the method must be either largerToSmaller or sMAPE.'
        if (len(np.shape(data))==2):
            data = self.traceToStack(data)
        end = np.shape(data)[-1]
        meanOdd = np.nanmean(data[:,:,0:end:2],2)
        meanEven = np.nanmean(data[:,:,1:end:2],2)

        maxMap = np.maximum(meanOdd, meanEven)
        minMap = np.minimum(meanOdd, meanEven)

        if (method == 'largerToSmaller'):
            alternans = maxMap / minMap
        elif (method == 'sMAPE'):
            alternans = abs(meanOdd-meanEven)/(meanOdd + meanEven)
        return alternans

    def getCV(self, activationMap, XY, varargin=None):
        # Analyse conduction velocity (CV) between pair (or pairs) of points.
        # By default, this returns the CV in pixels per frame, but can also
        # provide it in cm/s.
        #
        # IN:
        # activationMap - a single activation map.
        #
        # XY - a matrix of size n-by-4 encoding pairs of points between which CV is measured using the provided activation map.
        # Each row corresponds to an origin and target point (columns
        # are: rowFrom, columnFrom, rowTo, columnTo).
        #
        # varargin may optionally contain a two-numbers parameter allowing specification of spatial (how many mm is a single pixel side) and temporal resolution (in frames per second). If not given, the output is in pixels/frame, otherwise in cm/s.
        #
        # OUT:
        # cv -  a vector of conduction velocities, one per row of XY.
        XY = np.array(XY)
        nRows = np.shape(XY)[0]
        #CV = np.zeros(nRows)
        for iRow in range(0, nRows):
            fromRow = XY[iRow,0]
            fromCol = XY[iRow,1]
            toRow = XY[iRow,2]
            toCol = XY[iRow,3]
            distDiff = np.sqrt((fromRow-toRow)**2+(fromCol-toCol)**2)
            timeDiff = activationMap[toRow-1,toCol-1] - activationMap[fromRow-1, fromCol-1]

        cv = distDiff/timeDiff

        # Potential recasling to cm/s
        if (varargin):
            assert len(varargin)==2, 'the scaling vector'
            distancePerPixel = varargin[0]
            fps = varargin[1]
            cv = (cv * distancePerPixel*fps)/10.0
        return cv


    def getLocalCV(self, activationMap, baylyNeighbourhood, varargin=None):
        # Performs local estimation of CV using Bayly's method (doi: 10.1109/10.668746).
        #
        # IN:
        # activationMap - activation map from which a vector field is
        # obtained.
        #
        # baylyNeighbourhood - distance around a point that is
        # considered when fitting the Bayly polynomial.
        #
        # varargin{1} - if defined, gives maximum length of an arrow
        # (the longer ones are discarded).
        #
        # varargin{2} - if defined, contains the index of figure in
        # which the CV field is drawn. If not defined, no figure is
        # produced.
        #
        # varargin{3} - if defined, the output path of storage of
        # varargin{1}.
        #
        # OUT:
        # xyuv - a n-by-4 matrix, where n is number of pixels and columns correspond
        # to x,y,u,v:  x,y gives indices of row and column, with u,v
        # corresponding to dx,dy. Mind that this is in row/column
        # coordinates - if plotting via quiver (in standard x-y
        # coordinates), this needs to be altered slightly, see the code
        # at the end of the function.
        maxDistance = np.inf
        if (len(varargin)>=1):
            maxDistance = varargin[0] # maximum arrow length that is not discarded
        figureNumber = ''
        if (len(varargin)>=2):
            figureNumber = varargin[1] # figure number
        foutName = ''
        if (len(varargin)>=3):
            foutName = varargin[2] # output PNG file name

        ## Bayly CV estimation
        nRows = np.shape(activationMap)[0]
        nCols = np.shape(activationMap)[1]
        Z = np.zeros((nRows*nCols, 2)).astype(int)
        count = 0
        for i in range(0, nRows):
            for j in range(0, nCols):
                Z[count, 0] = j+1
                Z[count, 1] = i+1
                count = count + 1
        activationTimes = activationMap.reshape((nRows*nCols), order='F')
        xyuv = np.zeros((len(Z),4))

        for iPoint in range(0, len(Z)):
            # for each point, find points nearby.
            thisPoint = Z[iPoint,:]
            distances = cdist([thisPoint], Z)
            whereNeighbours = np.where((distances>=0)&(distances <=baylyNeighbourhood))[1]
            locationsNeighbours = Z[whereNeighbours,:]
            neighbourActivationTimes = activationTimes[whereNeighbours]
            sf = self.polyfit22(locationsNeighbours[:,0],locationsNeighbours[:,1],neighbourActivationTimes)
            x = np.arange(min(locationsNeighbours[:,0]), max(locationsNeighbours[:,0])+1)
            y = np.arange(min(locationsNeighbours[:,1]), max(locationsNeighbours[:,1])+1)
            coeffs = sf[0]
            x = thisPoint[0]
            y = thisPoint[1]
            dx = coeffs[1] + 2*coeffs[3]*x + coeffs[4]*y
            dy = coeffs[2] + 2*coeffs[5]*y + coeffs[4]*x
            xyuv[iPoint, :] = [x, y, dx/(dx*dx+dy*dy), dy/(dx*dx + dy*dy)]

        # We get rid of dx, dy, which are too long
        arrowLengths = np.sqrt(xyuv[:,2]*xyuv[:,2] + xyuv[:,3]*xyuv[:,3])
        idx = np.where(arrowLengths>maxDistance)
        xyuv[idx,2] = np.NaN
        xyuv[idx,3] = np.NaN

        ## Optional plotting of the quiver. Given that axes are
        # different between x/y (and for Matlab row/columns, [1,1] is
        # top left, while for common set of axes it is bottom left, so
        # we do complement for the 2nd parameter in quiver and we flip
        # the sign of the fourth parameter in quiver)

        if (figureNumber):
            plt.figure(figureNumber)
            plt.quiver(xyuv[:,1], 16-xyuv[:,0], xyuv[:,3], -xyuv[:,2])
            plt.axis([-5, 20, -5, 20])
            if (foutName):
                plt.savefig(foutName)
            else:
                plt.show()


    def plotActivationMap(self,activationMap, varargin=None):
        # This function plots a contour map of activation.
        #
        # IN:
        # activationMap - a single activation map from which contours are
        # obtained.
        #
        # varagin - optionally, the index of figure used for this purpose may be given - if not
        # specified, a new figure is opened and used.
        figureNumber = []
        levelGranularity = 1

        if (varargin):

            if type(varargin) is list:
                figureNumber = varargin[0]
                if (len(varargin) >= 2):
                    levelGranularity = varargin[1]
            else:
                figureNumber = varargin
                assert (figureNumber > 0), 'the second parameter (figure number) must be a positive integer'

        if not (figureNumber):
            plt.figure()
        else:
            plt.figure(figureNumber)
        end = np.round(np.nanmax(np.flipud(activationMap - np.nanmin(np.nanmin(activationMap)))))
        level = np.arange(0,end,levelGranularity).astype(int)

        plt.contourf(np.flipud(activationMap - np.nanmin(np.nanmin(activationMap))), levels=level, extend='both')
        cs = plt.contour(np.flipud(activationMap - np.nanmin(np.nanmin(activationMap))),levels=level,extend='both',linewidths=0.5,colors='k')
        plt.clabel(cs, levels=level, fmt='%d', colors='k')

        plt.clabel(cs, fmt='%d')
        if not os.path.exists('imOut'):
            os.mkdir('imOut')
        plt.savefig('imOut/sampleActivation.png')

##### HELPER FUNCTIONS #####

    def applyMask(self, imStackIn, mask):
        # This function may be used to apply a binary mask to each frame of a 3D stack, setting pixels-to-be-discarded to NaN.
        # This may be useful to get rid of empty space around the image
        # of the heart, space around a Petri dish for cultures, etc.
        #
        # IN:
        # imStackIn - a 3D stack representing a recording.
        #
        # mask - a binary mask of the same size as a single frame of
        # imStackIn. Where mask==0, the corresponding pixels are set to
        # NaN.
        #
        # OUT:
        # imStackOut - imStackIn where zero elements in mask are set to
        # NaN for each frame.

        imStackOut = imStackIn
        for iFrame in range(0, np.shape(imStackIn)[2]):
            frame = imStackIn[:,:,iFrame]
            frame[mask == False] = np.nan
            imStackOut[:,:,iFrame] = frame
        return imStackOut

    def binningVec(self,img,logfactor):
        # A function which carries out spatial binning for an image
        # (replacing each bin with the average of the pixels found in
        # it).
        #
        # IN:
        # img - an image to be spatially binned.
        #
        # logfactor - base-two logarithm of the binning factor - this must be a
        # positive integer (i.e., logfactor of 1 leads to 2-by-2
        # binning, 2 to 4-by4, 3 to 8-by-8, etc.). Note that while this
        # helper function requires a logarithm of the binning factor,
        # the data-reading functions readFolder and readTifStack take
        # the binning factor directly, computing its logarithm
        # internally, so that the user doesn't have to take care of
        # that.
        #
        # OUT:
        # binnedImage - the source image after spatial binning.
        binnedImage = img.astype('float')
        for j in range(0, int(logfactor)):
            end = len(binnedImage[:,0])
            oddRows = binnedImage[0:end:2]
            evenRows = binnedImage[1:end:2,:]
            rowBinned = (oddRows + evenRows) / 2.0
            end = len(rowBinned[0,:])
            oddCols = rowBinned[:,0:end:2]
            evenCols = rowBinned[:,1:end:2]
            binnedImage = self.roundHalfUpInt(((oddCols + evenCols)/2.0))
        return binnedImage

    def combGetMinima(self, signalTrace, bcl, varargin=None):
        # A function searching for minima in traces from cardiac preparations with known activation pattern.
        # The function can be naturally also used to extract signal maxima when
        # the source trace is inverted.
        #
        # IN:
        # signalTrace - a vector containing signal, such as calcium transients or action potentials
        # (e.g., intensity of a pixel in optical mapping, or an electrophysiological recording).
        #
        # bcl - the basic cycle length (number of frames between two activations).
        #
        # varargin{1} - The first optional parameter is the refinementWidth parameter for comb algorithm
        # (in ms/samples - the radius of local search around comb
        # teeth). Default is 10.
        #
        # varargin{2} - The second optional parameter is the custom comb that is used instead of the
        # regularly placed one. It should be an increasing vector of
        # numbers, where the first element determines the last possible
        # position of the first minimum to be searched for, and the
        # subsequent elements give further indices of minima. E.g.
        # using [30 180, 330, 400] means that the algorithm will search
        # for 3 minima that are 150 frames apart, and one that is
        # further 100 frames after the last previous one, and the first
        # minimum is to be placed between frame 1 and frame 30. See
        # sampleScript_5_S1S2 for an example of custom comb. If this parameter is defined, then
        # the parameter ‘bcl’ is not used within the code (it still has to be provided,
        # but it can be any number).
        #
        #
        # OUT:
        # vectorMinima - the vector of local minima in the signal, which are
        # approximately bcl ms (or samples) apart.

        ## Default parameter initialisation and processing of extra inputs.
        refinementWidth = 10
        customComb = []

        if (len(varargin) == 1):
            if varargin[0] is not None:
                refinementWidth = varargin[0]
        elif (len(varargin) == 2):
            if varargin[0] != []:
                refinementWidth = varargin[0]
            if varargin[1] != []:
                customComb = varargin[1]
        ## Comb positioning
        nFrames = len(signalTrace)

        if not customComb: # if no custom comb is given, we use a regular one based on bcl.
            candidateMinimaMean = np.full(bcl,np.inf)
            for iStart in range(0, bcl):
                candidateMinima = np.arange(iStart,nFrames,bcl,dtype=int)
                candidateMinimaMean[iStart] = np.average(signalTrace[candidateMinima])
        else: # when custom comb is provided, we slide it up to the first element
            candidateMinimaMean = np.full(customComb[0],np.inf)
            combNoOffset = customComb - customComb[0] # comb with first element subtracted so it starts at 0
            for iStart in range(0, customComb[0]):
                candidateMinima = iStart + combNoOffset
                candidateMinimaMean[iStart] = np.average(signalTrace[candidateMinima])

        # The start leading to minimal average value is used as the indicator of minima.
        bestMinima = np.where(candidateMinimaMean==candidateMinimaMean.min())[0][-1]
        if not customComb:
            vectorMinima = np.arange(bestMinima,nFrames,bcl)
        else:
            vectorMinima = bestMinima + combNoOffset

        ## Comb refinement
        # The vector of teeth (~near-minima) is refined by looking refinementWidth to the
        # left and right, taking the actual minimum there.
        vectorMinimaRefined = np.zeros(np.shape(vectorMinima))
        for iMinimum in range(0,len(vectorMinima)):
            leftBound = np.max([vectorMinima[iMinimum]-refinementWidth,0])
            rightBound = np.min([vectorMinima[iMinimum] + refinementWidth,len(signalTrace)])
            x= signalTrace[leftBound:rightBound+1]
            actualMinimum = np.where(x ==x.min())[0][-1]
            vectorMinimaRefined[iMinimum] = leftBound + actualMinimum
        vectorMinima = vectorMinimaRefined
        return vectorMinima

    def getAPD(self, time, activationSignal, level, varargin=None):
        # The function returns the duration of an action potential (or calcium transient)
        # at the desired level of repolarization.
        # The function assumes the action potential/calcium transient to be upward-pointing (i.e. peak activation
        # is more positive than resting value) - if
        # used on downward-pointing signal, this needs to be inverted first.
        #
        # IN:
        # time - the time vector for the activation trace.
        #
        # activationSignal - the trace of a single action potential or
        # calcium transient (or any similar signal).
        #
        # level - the level of recovery at which the duration is to be
        # measured. This is scaled between 0 and 1 (i.e., for APD80,
        # the value 0.8 is to be used).
        #
        # varargin{1} - a string encoding the method of baseline
        # estimation; see the documentation of SCUMS.analyseRegularPacing,
        # the 'parameters.baselineDefinition' parameter.
        #
        # varargin{2} - a string encoding the method of object
        # detection/selection; see the documentation of
        # SCUMS.analyseRegularPacing, the 'parameters.objectDetection'
        # parameter.
        #
        # varargin{3} - when 'augmented' objectDetection is used, this
        # parameter is used to pass the time of peak upstroke time
        # within the single activation provided in 'activationSignal'.
        # When more objects are above the threshold determined by
        # 'level', the one closest to this parameter is used.
        #
        # OUT:
        #
        # apd - the duration of the action potential/calcium transient
        # at the given level of recovery. The function uses interpolation to get sub-frame resolution.
        #
        # timeRecovery - the time of the end of the action
        # potential/calcium transient at the given repolarization
        # level. The function uses interpolation to get sub-frame resolution.
        activationSignal = activationSignal.astype('float')
        baselineDefinition = 'first'
        objectDetection = 'first'

        if (len(varargin)>=1):
            baselineDefinition = varargin[0]
        if (len(varargin)==2):
            objectDetection = varargin[1]
        if (len(varargin)==3):
            peakDiffActivation = varargin[2]

        baseline = activationSignal[0]
        if (baselineDefinition == 'firstlast'):
            baseline = np.average([activationSignal[0],activationSignal[-1]])

        thresh = baseline + (1-level)*(np.max(activationSignal)-baseline)
        s = measure.label(activationSignal > thresh, return_num=True)
        label = s[0]
        numObjects = s[1]
        pixelIdxList = {}
        for k in range(0, numObjects):
            pixelIdxList[k] = [i for i,v in enumerate(label) if v == (k+1)]

        # Now, the longest object is found among those above threshold.
        # It's almost always the first one, but if there is a lot of
        # noise (or pepper noise), it's not necessarily the case.
        ## We switch between taking the largest and first found object
        if (numObjects > 1) & (objectDetection == 'largest'):
            lengths = []
            for i in range(0, len(pixelIdxList)):
                lengths.append(len(pixelIdxList[i]))
            lengths = np.array(lengths)
            whereMax = np.where(lengths==max(lengths))[0][-1]
            interval = pixelIdxList[whereMax]
        elif (numObjects > 1) & (objectDetection == 'augmented'):
            assert peakDiffActivation, 'when augmented object detection used, peakDiffActivation parameter (3rd varargin) must be specified'
            distances = np.zeros((numObjects,1))
            for iObject in range(0, numObjects):
                pix = pixelIdxList[iObject]
                distances[iObject] = np.min(np.abs(pix-peakDiffActivation))
            whereMin = np.where(distances == distances.min())[-1]
            interval = pixelIdxList[whereMin]
        else:
            interval = pixelIdxList[0]

        # Linear interpolation of times when the signal crosses
        # threshold at the activation and recovery is carried out
        tStartInterval = time[interval[0]]
        tEndInterval = time[interval[-1]]

        ## Linear interpolation of the dt after the segment above threshold ends.
        tEndInterval = time[interval[-1]]
        try:
            tBeforeStartInterval = time[interval[0]-1]
            tAfterEndInterval = time[interval[-1]+1]
        except:
            apd = np.nan
            timeRecovery = np.nan
            return apd, timeRecovery

        # Processing interval before activation start
        vStartInterval = activationSignal[interval[0]]
        vBeforeStartInterval = activationSignal[interval[0]-1]
        howFarThreshBeforeStart = (vStartInterval-thresh)/(vStartInterval-vBeforeStartInterval)
        interpAddingStart = (tStartInterval-tBeforeStartInterval) * howFarThreshBeforeStart
        # Processing interval after recovery
        vEndInterval = activationSignal[interval[-1]]
        vAfterEndInterval = float(activationSignal[np.max(interval)+1])
        howFarThreshAfterEnd = (thresh - vEndInterval) / (vAfterEndInterval- vEndInterval)
        interpAddingEnd = (tAfterEndInterval - tEndInterval) * howFarThreshAfterEnd
        apd = time[interval[-1]] - time[interval[0]] + interpAddingStart + interpAddingEnd

        timeRecovery = time[interval[-1]] + interpAddingEnd

        return apd, timeRecovery


    def getHalfActivationTime(self, activationSignal, varargin=None):
        # The function returns the time at "half-activation" of an action potential
        # or calcium transient; taken with regards to amplitude (i.e. it corresponds
        # to the time when APD50 would be measured from when applied to action potentials).
        # It uses linear interpolation to get this half-activation time
        # "accurately" (as in not just the previous/subsequent frame).
        # The function assumes the action potential/calcium transient to be upward-pointing
        # (i.e. peak activation is more positive than resting value) - if used
        # on downward-pointing signal, this needs to be inverted first.
        #
        # IN:
        # activationSignal - the trace of a single action potential or
        # calcium transient (or any similar signal).
        #
        # varargin{1} - a string encoding the method of baseline
        # estimation; see the documentation of SCUMS.analyseRegularPacing,
        # the 'parameters.baselineDefinition' parameter.
        #
        # varargin{2} - a string encoding the method of object
        # detection/selection; see the documentation of
        # SCUMS.analyseRegularPacing, the 'parameters.objectDetection'
        # parameter.
        #
        # varargin{3} - when 'augmented' objectDetection is used, this
        # parameter is used to pass the time of peak upstroke time
        # within the single activation provided in 'activationSignal'.
        # When more objects are above the threshold determined by
        # 'level', the one closest to this parameter is used.
        #
        # OUT:
        # halfActivationTime - The time of half-activation of the signal
        # (i.e. time when 50# of the signal amplitude is achieved).
        # The function uses interpolation to get sub-frame resolution.
        baselineDefinition = 'first'
        objectDetection = 'first'
        peakDiffActivation = []
        if (varargin):
            baselineDefinition = varargin[0]
            if (len(varargin) >= 2):
                objectDetection = varargin[1]
            if (len(varargin) == 3):
                peakDiffActivation = varargin[2]

        try:
            baseline = activationSignal[0]
            if (baselineDefinition == 'firstlast'):
               baseline = np.average([activationSignal[0], activationSignal[-1]])

            peak = np.min(activationSignal)
            threshold = baseline - 0.5 * (baseline - peak)
            s = measure.label(activationSignal < threshold, return_num=True)
            label = s[0]
            numObjects = s[1]
            pixelIdxList = {}
            for k in range(0, numObjects):
                pixelIdxList[k] = [i for i,v in enumerate(label) if v == (k+1)]
            if (numObjects > 1) & (objectDetection == 'largest'):
                lengths = []
                for i in range(0, len(pixelIdxList)):
                    lengths.append(len(pixelIdxList[i]))
                lengths = np.array(lengths)
                whereMax = np.where(lengths==max(lengths))[0][-1]
                apFrames = pixelIdxList[whereMax]
            elif (numObjects > 1) & (objectDetection == 'augmented'):
                assert peakDiffActivation, 'when augmented object detection used, peakDiffActivation parameter (3rd varargin) must be specified'
                distances = np.zeros((numObjects,1))
                for iObject in range(0, numObjects):
                    pix = pixelIdxList[iObject]
                    distances[iObject] = np.min(np.abs(pix - peakDiffActivation))
                whereMin = np.where(distances == distances.min())[0][-1]
                apFrames = pixelIdxList[whereMin]
            else:
                apFrames = pixelIdxList[0]
            apFrames = np.array(apFrames).astype(int)
            afterHalfActivationTime = apFrames[0]
            preHalfActivationTime = apFrames[0]-1
            traceAfter = float(activationSignal[afterHalfActivationTime])
            tracePre = float(activationSignal[preHalfActivationTime])
            halfActivationTime = preHalfActivationTime + ((threshold - tracePre)/ (traceAfter - tracePre)) * (afterHalfActivationTime - preHalfActivationTime) + 1 # the last parentheses are 1 ms at the moment - written like this just for future generality.
            return halfActivationTime
        except:
            print('issue in getting half activation')
            halfActivationTime = np.nan
            return halfActivationTime

    def getImageForMask(self, imStack):
        # This function gives a high-contrast average image of the stack (averaged over
        # time). This can be used to draw a mask (manually or automatically).
        #
        # IN:
        # imStack - A stack representing the recording. The input can be of type
        # uint8, uint16, or double (that one is expected to be scaled between
        # 0 and 1 - if not, please convert the stack to the uint types).
        #
        # OUT:
        # imOut - A maximum-contrast image of the time-average of the
        # stack.
        avgIm = np.average(imStack,2)
        avgIm = avgIm - np.min(avgIm) # subtracting baseline
        maxIm = np.max(avgIm)
        if (avgIm.dtype == 'float64'):
            imOut = avgIm * 1/ maxIm
        elif (avgIm.dtype == 'uint8'):
            imOut = avgIm * 255 / maxIm
        elif (avgIm.dtype == 'uint16'):
            imOut = avgIm * 65535/maxIm
        else:
            exit('Stack type must be double, uint8, uint16.')
        return imOut

    def processSingleActivation(self, activationTrace,spikesPointDown,baselineDefinition,durationLevel,objectDetection,recordingClock,timeStart,peakDiffActivation,verbose):
        # Processes a single activation (action potential or calcium
        # transient), returning its properties/features - this is
        # mainly a helper function for SCUMS.analyseRegularPacing or
        # SCUMS.analyseSinglePass, and the parameters it uses are the
        # same as there.

        # Sorting where data will be stored within the
        # coordinates given by recordingClock.
        if (spikesPointDown):
            wherePeak = np.where(activationTrace==activationTrace.min())[0][0]
        else:
            wherePeak =np.where(activationTrace==activationTrace.max())[0][0]
        # Try-catch, because when signal is poor,
        # getHalfActivation often fails, and in that case,
        # we don't really want the other values either
        try:
            # Baseline - first index
            bas = activationTrace[0]
            if baselineDefinition == 'firstlast':
                bas = np.average([activationTrace[0], activationTrace[-1]])
            saBaseline = bas

            # Amplitude
            saAmplitude = abs(activationTrace[wherePeak] - bas)

            # Duration
            if (spikesPointDown):
                saDuration, saRecovery = self.getAPD(np.arange(0,len(activationTrace),1),-1*activationTrace, durationLevel, varargin=[baselineDefinition, objectDetection, peakDiffActivation])
            else:
                saDuration, saRecovery = self.getAPD(np.arange(0, len(activationTrace),1), activationTrace, durationLevel, varargin=[baselineDefinition, objectDetection, peakDiffActivation])
            # Activation - first extracted in local coordinates, converted to global
            # coordinates, and then it is taken relatively
            # to the clock
            if (spikesPointDown):
                activationLocal = self.getHalfActivationTime(activationTrace, [baselineDefinition, objectDetection, peakDiffActivation])
            else:
                activationLocal = self.getHalfActivationTime(-1*activationTrace, [baselineDefinition, objectDetection, peakDiffActivation])
            activationGlobal = timeStart + activationLocal
            whereLastSmaller = np.where(recordingClock<activationGlobal)
            if (not whereLastSmaller):
                exit('There is a problem with data, possibly with the rate of activation, or there may be an arrhythmia')
            whereLastSmaller = whereLastSmaller[0][-1]

            saActivation = activationGlobal-recordingClock[whereLastSmaller]

            if (not np.isnan(saRecovery)):
                recoveryGlobal = timeStart + saRecovery
                saRecovery = recoveryGlobal - recordingClock[whereLastSmaller]
            else:
                saRecovery = np.nan
        except:
            if (verbose):
                print('problem processing a trace of a pixel')
            saBaseline = np.nan
            saAmplitude = np.nan
            saDuration = np.nan
            saActivation = np.nan
            saRecovery = np.nan
        return saBaseline, saAmplitude, saDuration, saActivation, saRecovery


    def smoothTrace(self, tracePixel, filterWidth, varargin=None):
        # This function smooths (~denoises) a trace using Savitzky-Golay filter of
        # a given width. It may also remove drift/trend in data using
        # polynomial fitting.
        #
        # IN:
        # tracePixel - the trace to be smoothed.
        #
        # filterWidth - the width of the Savitzky-Golay filter (should
        # be an odd number - if an even one is provided, one is added
        # to it).
        #
        # varargin{1} - the optional parameter may contain the order of the polynomial
        # used to detrend the signal (useful for removing drift etc.).
        # Unlike standard methods which subtract the polynomial making
        # the signal approximately zero-centered on the y-axis, here we
        # re-add the average of the original trace to the zero-centered
        # trace, so that the information on signal baseline is not
        # lost.
        #
        # OUT:
        # traceSmoothed - tracePixel after smoothing.
        #
        # traceSmoothedDetrended - tracePixel after smoothing and
        # baseline subtraction.
        tracePixel = tracePixel[:]
        if varargin:
            detrendOrder = varargin
            assert (detrendOrder >0), 'Order of baseline subtraction must be larger than 0'
            # Below we fit a polynomial to the signal, which acts as a baseline.
            # Mu shenanigans serve to improve the numerical stability (see documentation of polyfit for mu)
            p,v = np.polyfit(np.arange(1,len(tracePixel)+1), tracePixel, detrendOrder, cov=True, full=False)
            polyBaseline = np.polyval(p,np.arange(1,len(tracePixel)+1))
            traceDetrend = tracePixel - polyBaseline + np.average(tracePixel)
            #traceDetrend = signal.detrend(tracePixel, type='linear')
        if (np.mod(filterWidth,2) == 0):
            filterWidth = filterWidth + 1
        if (filterWidth > 4): # if filter width is <= 4, it's not possible to carry out filtering. In this way, we can also forbid filtering (by setting smoothingParameter to 1)
            traceSmoothed = signal.savgol_filter(np.ravel(tracePixel),filterWidth,4)
            traceSmoothedDetrended = traceSmoothed # just so that this is defined as an output parameter if no detrending is used.
            if varargin:
                traceSmoothedDetrended = signal.savgol_filter(np.ravel(traceDetrend),filterWidth,4)
                return traceSmoothed, traceSmoothedDetrended
            else:
                return traceSmoothed
        else:
            traceSmoothed = tracePixel
            traceSmoothedDetrended = traceSmoothed
            if varargin:
                traceSmoothedDetrended = traceDetrend
                return traceSmoothed, traceSmoothedDetrended
            else:
                return traceSmoothed


    def traceToStack(self, signalTrace):
        # A minor helper function which converts a vector (either n-by-1 or 1-by-n)
        # to a “fake” 3D stack of 1-by-1-by-n – in this way, it is the same format
        # as a stack representing a video (albeit with 1x1 resolution) and it can
        # be fed to analyseRegularPacing. One can therefore use the functionality
        # of SCUMS with regards to baseline, amplitude, or duration features even
        # on single traces (coming e.g. from electrophysiological recordings).
        #
        # IN:
        # signalTrace - the trace to be converted into a "stack".
        #
        # OUT:
        # imStack - the resulting stack.
        imStack = np.zeros([1,1,len(signalTrace)])
        imStack[0,0,:] = signalTrace[:,0]
        return imStack

##### TECHNICAL HELPER #####
    def roundHalfUpInt(self, n):
        if (hasattr(n, '__len__')):
            if (len(np.shape(n)) == 1):
                for i in range(0, len(n)):
                    n[i] = math.floor((1*n[i] + 0.5)/1)
            elif (len(np.shape(n)) == 2):
                for i in range(0, len(n[:,0])):
                    for j in range(0, len(n[0,:])):
                        n[i,j] = math.floor((n[i,j]+0.5))
            return n.astype('uint16')
        else:
            return int(math.floor((1*n + 0.5)/1))

    def nanmean(self, a, axis=None, dtype=None, out=None, keepdims=np._NoValue):
        """
        Compute the arithmetic mean along the specified axis, ignoring NaNs.
        Returns the average of the array elements.  The average is taken over
        the flattened array by default, otherwise over the specified axis.
        `float64` intermediate and return values are used for integer inputs.
        For all-NaN slices, NaN is returned and a `RuntimeWarning` is raised.
        .. versionadded:: 1.8.0
        Parameters
        ----------
        a : arrayLike
            Array containing numbers whose mean is desired. If `a` is not an
            array, a conversion is attempted.
        axis : {int, tuple of int, None}, optional
            Axis or axes along which the means are computed. The default is to compute
            the mean of the flattened array.
        dtype : data-type, optional
            Type to use in computing the mean.  For integer inputs, the default
            is `float64`; for inexact inputs, it is the same as the input
            dtype.
        out : ndarray, optional
            Alternate output array in which to place the result.  The default
            is ``None``; if provided, it must have the same shape as the
            expected output, but the type will be cast if necessary. See
            `ufuncs-output-type` for more details.
        keepdims : bool, optional
            If this is set to True, the axes which are reduced are left
            in the result as dimensions with size one. With this option,
            the result will broadcast correctly against the original `a`.
            If the value is anything but the default, then
            `keepdims` will be passed through to the `mean` or `sum` methods
            of sub-classes of `ndarray`.  If the sub-classes methods
            does not implement `keepdims` any exceptions will be raised.
        Returns
        -------
        m : ndarray, see dtype parameter above
            If `out=None`, returns a new array containing the mean values,
            otherwise a reference to the output array is returned. Nan is
            returned for slices that contain only NaNs.
        See Also
        --------
        average : Weighted average
        mean : Arithmetic mean taken while not ignoring NaNs
        var, nanvar
        Notes
        -----
        The arithmetic mean is the sum of the non-NaN elements along the axis
        divided by the number of non-NaN elements.
        Note that for floating-point input, the mean is computed using the same
        precision the input has.  Depending on the input data, this can cause
        the results to be inaccurate, especially for `float32`.  Specifying a
        higher-precision accumulator using the `dtype` keyword can alleviate
        this issue.
        Examples
        --------
        >>> a = np.array([[1, np.nan], [3, 4]])
        >>> self.nanmean(a)
        2.6666666666666665
        >>> self.nanmean(a, axis=0)
        array([2.,  4.])
        >>> self.nanmean(a, axis=1)
        array([1.,  3.5]) # may vary
        """
        arr, mask = self.ReplaceNan(a, 0)
        if mask is None:
            return np.mean(arr, axis=axis, dtype=dtype, out=out, keepdims=keepdims)
        if dtype is not None:
            dtype = np.dtype(dtype)
        if dtype is not None and not issubclass(dtype.type, np.inexact):
            raise TypeError("If a is inexact, then dtype must be inexact")
        if out is not None and not issubclass(out.dtype.type, np.inexact):
            raise TypeError("If a is inexact, then out must be inexact")
        cnt = np.sum(~mask, axis=axis, dtype=np.intp, keepdims=keepdims)
        tot = np.sum(arr, axis=axis, dtype=dtype, out=out, keepdims=keepdims)
        avg = self.DivideByCount(tot, cnt, out=out)
        isbad = (cnt == 0)
        return avg
    def ReplaceNan(self, a, val):
        """
        If `a` is of inexact type, make a copy of `a`, replace NaNs with
        the `val` value, and return the copy together with a boolean mask
        marking the locations where NaNs were present. If `a` is not of
        inexact type, do nothing and return `a` together with a mask of None.
        Note that scalars will end up as array scalars, which is important
        for using the result as the value of the out argument in some
        operations.
        Parameters
        ----------
        a : array-like
            Input array.
        val : float
            NaN values are set to val before doing the operation.
        Returns
        -------
        y : ndarray
            If `a` is of inexact type, return a copy of `a` with the NaNs
            replaced by the fill value, otherwise return `a`.
        mask: {bool, None}
            If `a` is of inexact type, return a boolean mask marking locations of
            NaNs, otherwise return None.
        """
        a = np.asanyarray(a)
        if a.dtype == np.object_:
            # object arrays do not support `isnan` (gh-9009), so make a guess
            mask = np.notEqual(a, a, dtype=bool)
        elif issubclass(a.dtype.type, np.inexact):
            mask = np.isnan(a)
        else:
            mask = None
        if mask is not None:
            a = np.array(a, subok=True, copy=True)
            np.copyto(a, val, where=mask)
        return a, mask
    def DivideByCount(self, a, b, out=None):
        """
        Compute a/b ignoring invalid results. If `a` is an array the division
        is done in place. If `a` is a scalar, then its type is preserved in the
        output. If out is None, then then a is used instead so that the
        division is in place. Note that this is only called with `a` an inexact
        type.
        Parameters
        ----------
        a : {ndarray, numpy scalar}
            Numerator. Expected to be of inexact type but not checked.
        b : {ndarray, numpy scalar}
            Denominator.
        out : ndarray, optional
            Alternate output array in which to place the result.  The default
            is ``None``; if provided, it must have the same shape as the
            expected output, but the type will be cast if necessary.
        Returns
        -------
        ret : {ndarray, numpy scalar}
            The return value is a/b. If `a` was an ndarray the division is done
            in place. If `a` is a numpy scalar, the division preserves its type.
        """
        with np.errstate(invalid='ignore', divide='ignore'):
            if isinstance(a, np.ndarray):
                if out is None:
                    return np.divide(a, b, out=a, casting='unsafe')
                else:
                    return np.divide(a, b, out=out, casting='unsafe')
            else:
                if out is None:
                    return a.dtype.type(a / b)
                else:
                    # This is questionable, but currently a numpy scalar can
                    # be output to a zero dimensional array.
                    return np.divide(a, b, out=out, casting='unsafe')

    def polyfit22(self, x, y, z):
        coeffs = np.ones(6)
        a = np.zeros((coeffs.size, x.size))
        a[0] = (coeffs[0] * x**0 * y**0).ravel()
        a[1] = (coeffs[1] * x**1 * y**0).ravel()
        a[2] = (coeffs[2] * x**0 * y**1).ravel()
        a[3] = (coeffs[3] * x**2 * y**0).ravel()
        a[4] = (coeffs[4] * x**1 * y**1).ravel()
        a[5] = (coeffs[5] * x**0 * y**2).ravel()
        # do leastsq fitting and return leastsq result
        return np.linalg.lstsq(a.T, np.ravel(z), rcond=None)

    def polyval22(self, coeffs, x, y, xydata, zdata):

        X, Y = np.meshgrid(x,y)
        Z = coeffs[0] * X**0 * Y**0 + coeffs[1] * X**1 * Y**0 + \
            coeffs[2] * X**0 * Y**1 + coeffs[3] * X**2 * Y**0 + \
            coeffs[4] * X**1 * Y**1 + coeffs[5] * X**0 * Y**2
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(X, Y, Z,linewidth=0, antialiased=False)
        xdata = xydata[:,0]
        ydata = xydata[:,1]
        ax.scatter3D(xdata, ydata, zdata);
        plt.show()

    def polyfit2d_v2(self, x, y, z, order=2):
        ncols = (order + 1)**2
        G = np.zeros((x.size, ncols))
        # ij = itertools.product(range(order+1), range(order+1))
        counter = 0
        for i in range(0, order+1):
            for j in range(0, order+1):
                G[:, counter] = x**i * y**j
                counter = counter + 1
        # for k, (i,j) in enumerate(ij):
        #     G[:,k] = x**i * y**j
        output = np.linalg.lstsq(np.ravel(G).T, np.ravel(z))
        m = output[0]
        return m
