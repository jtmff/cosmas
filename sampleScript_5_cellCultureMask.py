import os
import numpy as np
import matplotlib.pyplot as plt
from COSMAS import COSMAS
from imageio import imwrite, imread
from scipy.io import loadmat
cosmas = COSMAS()

# This script demonstrates how to a) load a stack of data stored in a .mat
# file, rather than in image files, b) apply binary mask to it to discard
# image background. In this case, this is applied to remove the
# surroundings of a Petri dish with a cell culture. The script ultimately
# plots the conduction map (based on calcium mapping data).

class paramReading:
    extension = ''
    binningFactor = ''
class parameters:
    spikesPointDown = True
    baselineDefinition = ''

data = loadmat('data/cultureStack/stack.mat')
imStackCulture = data['imStackCulture'].astype('float')

# This is what the average frame looks like - it is stored and then its
# copy can be painted over with black color in the areas that are not to be
# processed.
imForMaskCulture = cosmas.getImageForMask(imStackCulture.astype("uint16"))
imwrite('masks/cultureAvg_python.png',imForMaskCulture.astype('uint8'))

# Now we applied red color to parts of the image we want to ignore
# (assuming the rest is grayscale) - binary mask is then obtained by
# keeping only pixels that were not painted in red color.
# In this case, the mask was drawn manually, but in general, it can be
# often also constructed automatically based on imForMaskCulture (e.g. via
# intensity thresholding)
mask = imread('masks/cultureMask.png')
maskBinary = (mask[:,:,0] == mask[:,:,1]) # If red==green (i.e. no red painting), we keep image.

imStackCultureMasked = cosmas.applyMask(imStackCulture, maskBinary)

parameters.spikesPointDown = False
parameters.baselineDefinition = 'first'
parameters.smoothingParameter = 21
baseline, amplitude, duration, activationMaps, recoveryMaps, clockFiltered = cosmas.analyseRegularPacing(imStackCultureMasked, 200, parameters) # 1 Hz, but 200 fps recording
am = activationMaps.mapMean
cosmas.plotActivationMap(am, [1,1,True])

plt.show()
