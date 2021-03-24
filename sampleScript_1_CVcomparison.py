import os
import numpy as np
import matplotlib.pyplot as plt
from COSMAS import COSMAS
cosmas = COSMAS()

# This script takes two recordings of voltage mapping in the same heart,
# before and after application of isoproterenol, showing an increase in
# conduction velocity under isoproterenol.

# Setting up classes for input parameters
class paramReading:
    extension = ''
    binningFactor = ''
class parameters:
    spikesPointDown = True
    baselineDefinition = ''
    objectDetection = ''
# Reading png images in the 'controlV' folder, performing 2x2 binning
paramReading.extension = '.png'
paramReading.binningFactor = 2
imStackV = cosmas.readFolder('data/controlV', paramReading)[0]

# Setting parameters of processing and processing itself
parameters.spikesPointDown = True
parameters.baselineDefinition = 'firstlast'
parameters.objectDetection = 'augmented'
bcl = 150 # a stimulus was applied every 150 ms
# Analysing the data
baseline, amplitude, duration, activationMaps, recoveryMaps, clockFiltered = cosmas.analyseRegularPacing(imStackV, bcl, parameters)

# Plotting the mean activation map
cosmas.plotActivationMap(activationMaps.mapMean)
if not os.path.exists('imOut'):
    os.mkdir('imOut')
plt.savefig('imOut/sampleActivation.png')

# Reading and analysing the folder with recording of the heart with isoproterenol
imStackViso, at = cosmas.readFolder('data/isoV', paramReading)
baselineISO, amplitudeISO, durationISO, activationMapsISO, recoveryMapsISO, clockFilteredISO = cosmas.analyseRegularPacing(imStackViso, bcl, parameters)

# Taking mean activation maps for control and iso
am = activationMaps.mapMean
amISO = activationMapsISO.mapMean

plt.figure(1)
plt.subplot(1,2,1)
end = np.round(np.max(np.flipud(am - np.min(np.min(am)))))
level = np.arange(0,end).astype(int)
#cs = plt.contour(np.flipud(am - np.min(np.min(am))),levels=level,extend='both')
plt.contourf(np.flipud(am - np.nanmin(np.nanmin(am))), levels=level, extend='both')
cs = plt.contour(np.flipud(am - np.nanmin(np.nanmin(am))),levels=level,extend='both',linewidths=0.5,colors='k')
plt.clabel(cs, levels=level, fmt='%d', colors='k')

plt.title('Control')
end = np.round(np.max(np.flipud(amISO - np.min(np.min(amISO)))))
level = np.arange(0,end).astype(int)
plt.subplot(1,2,2)
plt.contourf(np.flipud(amISO - np.nanmin(np.nanmin(amISO))), levels=level, extend='both')
cs = plt.contour(np.flipud(amISO - np.nanmin(np.nanmin(amISO))),levels=level,extend='both',linewidths=0.5,colors='k')
plt.clabel(cs, levels=level, fmt='%d', colors='k')
plt.title('ISO')

# Measuring conduction velocity between row=13,column=2 and row=1,column=2.
# The factor of 7/16 is the pixel side length in mm (FOV is 7-by-7mm and we
# have 16-by-16 pixels after the binning), and 1000 is the sampling
# frequency (1000 Hz).
cv1 = cosmas.getCV(am, [[13, 2, 1, 2]], [7/16, 1000])
cv2 = cosmas.getCV(amISO, [[13, 2, 1, 2]], [7/16, 1000])
print('Relative speedup: '+str(cv2/cv1))

plt.show()
