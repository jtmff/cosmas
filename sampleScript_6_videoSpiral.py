import cv2 # install via: pip3 install opencv-python
import numpy as np
from COSMAS import COSMAS
import matplotlib.pyplot as plt
from skimage import morphology
import imagesc
import code
cosmas = COSMAS()

# This scripts shows how to a) read data of a spiral wave mapping using
# off-axis illumination imaging from an mp4 video, b) plot the activation
# map and show how standard deviation of duration may be used to get rid of
# spurious signal areas.
class parameters:
    spikesPointDown = True
    baselineDefinition = ''
    objectDetection = ''
    baselineSubtractionOrder = ''
    smoothingParameter = ''

binningFactor = 4 # we do 4x4 binning on the unnecessarily high-res video.

v = cv2.VideoCapture('data/spiral_wave.mp4')
width = v.get(cv2.CAP_PROP_FRAME_WIDTH)
height = v.get(cv2.CAP_PROP_FRAME_HEIGHT)
num_frames = int(v.get(cv2.CAP_PROP_FRAME_COUNT))
stack = np.zeros((int(height/binningFactor), int(width/binningFactor), num_frames))

for i in range(0, num_frames):
    print(i)
    ret, frame = v.read()
    frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    frame = cosmas.binningVec(frame, np.log2(binningFactor))
    stack[:,:,i] = frame
parameters.baselineSubtractionOrder = 3
parameters.smoothingParameter = 21
parameters.spikesPointDown = False
parameters.baselineDefinition = 'firstlast'
parameters.objectDetection = 'first'
parameters.verbose = False

bcl = 34
msPerFrame = 10 # culture imaged at 100 fps

# Feature extraction is performed (given noisy signal and pixels with
# irregular activaton, there are some instances of 'issue in getting half
# activation' - this is just a warning, not a problematic error. The
# resulting pixels will have the feature values 'NaN'.
baseline, amplitude, duration, activationMaps, recoveryMaps, clockV = cosmas.analyseRegularPacing(stack, bcl, parameters)

# activation map of the spiral wave is plotted
cosmas.plotActivationMap(activationMaps.mapMean*msPerFrame,[1,20])

## Making a mask for visualisation cleaning based on postprocessing of duration maps
# (a manually drawn mask might be more accurate nevertheless, this is only
# a proof of concept)
#
# The standard deviation is computed over the 3rd dimension of the stack of
# signal durations, showing high values in areas of irregular activation,
# while giving low values where the main spiral is present.
mapStd = np.nanstd(duration.maps*msPerFrame, axis=2)

f, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.imshow(mapStd)
ax1.set_title('Std of duration (over wave passes) map')

maskToDiscard = ((np.isnan(mapStd)) | (mapStd >3*msPerFrame) ).astype(int)

ax2.imshow(maskToDiscard,cmap='Greys',  interpolation='nearest')
ax2.set_title('for removel - before filtering')

# Fill in the holes
img = np.pad(maskToDiscard, pad_width=1)
seed = np.ones_like(img)*255
img[ : ,0] = 0
img[ : ,-1] = 0
img[ 0 ,:] = 0
img[ -1 ,:] = 0
seed[ : ,0] = 0
seed[ : ,-1] = 0
seed[ 0 ,:] = 0
seed[ -1 ,:] = 0
fill = morphology.reconstruction(seed, img, method='erosion')
maskToDiscard = fill[1:-1, 1:-1]

# Perform the opening and closing operations
# morphological operations
str = morphology.ball(1)
strLarge = morphology.ball(3)
kernel = str[:,:,1]
kernelLarge = strLarge[:,:,1]
maskToDiscard = cv2.morphologyEx(maskToDiscard, cv2.MORPH_OPEN, kernel)
maskToDiscard = cv2.morphologyEx(maskToDiscard, cv2.MORPH_CLOSE, kernelLarge)
ax3.imshow(maskToDiscard, cmap='Greys',  interpolation='nearest')
ax3.set_title('for removel - after filtering')

# filteringthe activation map to set 'for removal' pixels to NaN
actMapFiltered = activationMaps.mapMean
idx = np.where(maskToDiscard == 1)
actMapFiltered[idx] = np.nan

# plotting activation map after postprocessing
cosmas.plotActivationMap(actMapFiltered*msPerFrame, [3, 20])

plt.show()
