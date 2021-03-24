import os
import numpy as np
import matplotlib.pyplot as plt
from COSMAS import COSMAS
cosmas = COSMAS()

# This script shows how to extract properties such as APD from single-trace
# recordings. While such a trace here was extracted as a spatial average of
# imaging experiment, electrophysiological recordings may be naturally
# processed in a similar fashion.
# Setting up classes for input parameters
class paramReading:
    extension = ''
    binningFactor = ''
class parameters:
    spikesPointDown = True
    baselineDefinition = ''

paramReading.extention = '.png'
paramReading.binningFactor = 2
# the second output parameter of readFolder is the signal trace produced as a
# result of spatially averaging the field of view in each frame.
imStackViso, at = cosmas.readFolder('data/isoV', paramReading)

# the trace (a vector) is converted to a "3D stack" with dimensions
# 1-by-1-by-n so that analyseRegularPacing may be directly applied
stackTrace = cosmas.traceToStack(at)
parameters.spikesPointDown = True # analysing membrane potential, with dye intensity reducing upon membrane depolarisation
parameters.baselineDefinition = 'firstlast'

# analysing the trace
baseline, amplitude, duration, activationMaps, recoveryMaps, clockFiltered = cosmas.analyseRegularPacing(stackTrace, 150, parameters)

# plotting the trace and showing the mean APD in the title
# (Please mind that taking average APD in a spatially averaged recording
# has only limited relevance and it may be often more representative to
# measure APD for each pixel and then average that - especially for preps
# with slow conduction such as certain cell cultures, the APD of
# whole-FOV-average is strongly modulated by conduction velocity. While a
# recording with infinite conduction velocity would have the same
# single-pixel and whole-FOV APD, if the conduction velocity is low, the
# temporally shifted action potentials combine into one that is possibly
# substantially prolonged compared to the single-pixel one).
plt.figure(1)
plt.plot(at)
plt.xlabel('Frame (no.)')
plt.ylabel('Voltage dye intensity')
title = 'APD75 observed in the recording: {:.4}'.format(duration.dataMean[0][0])
plt.title(title) # APD75 is the default level for APD extraction.
plt.show()
