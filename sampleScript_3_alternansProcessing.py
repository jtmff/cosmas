import os
import numpy as np
import matplotlib.pyplot as plt
from COSMAS import COSMAS
cosmas = COSMAS()

# This script shows how to analyse calcium transient (CaT) alternans in a rapidly-paced heart. It
# visualizes and quantifies the discordant alternans present.
class paramReading:
    extension = ''
    binningFactor = ''
class parameters:
    spikesPointDown = True
    baselineDefinition = ''

paramReading.extension = '.png'
paramReading.binningFactor = 2
imStackCaAlt, atCaAlt = cosmas.readFolder('data/alternansDiscordantCa_bcl85', paramReading)

parameters.spikesPointDown = False
parameters.baselineDefinition = 'first' #  if we want to measure alternans in amplitude of CaT upstroke, we measure baseline as the first element in each segmented calcium transient.
baselineCaAlt, amplitudeCaAlt, durationCaAlt, activationMapsCaAlt, recoveryMapsCaAlt, clockFilteredCaAlt = cosmas.analyseRegularPacing(imStackCaAlt, 85, parameters)

alternansMap = cosmas.getAlternans(amplitudeCaAlt.maps,'sMAPE')

# We plot the alternans map - cold colours correspond to nodal lines, where
# no alternans is present
fig = plt.figure(1)
ax= plt.gca()
im = ax.imshow(alternansMap)
ax.set_aspect('auto')
ax.set_title('alternans map, mean alternans = '+str(np.average(alternansMap[:])))
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel("", rotation=-90, va="bottom")

# As an additional visualization, we show the mean CaT for even/odd wave
# passes - where these differ, alternans is present
fig2 = plt.figure(2)
fig2.suptitle('CaT amplitude for even/odd waves')
ax = plt.subplot(1,2,1)
im = ax.imshow(amplitudeCaAlt.mapMeanEven)
ax.set_aspect('auto')
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel("", rotation=-90, va="bottom")
ax = plt.subplot(1,2,2)
ax.imshow(amplitudeCaAlt.mapMeanOdd)
ax.set_aspect('auto')
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel("", rotation=-90, va="bottom")

plt.show()

# Note: average alternans should be ideally taken as an average of an
# alternans map. An alternative and potentially appealing way is:
alternansData = cosmas.getAlternans(amplitudeCaAlt.data, 'sMAPE');
# , where alternans on average CaT per-wave-pass is measured. However, this
# is confounded by the discordant nature of alternans; just
# taking averages between even and odd averages for calcium transient is
# overshadowing the discordance. If unclear, imagine a discordant alternans
# where the left half of recording is in the opposite phase of the right
# half (and the sides are otherwise identical) - then there is no
# difference betwen beats in the calcium transient amplitude when averaged
# over the field of view, yet there can be pronounced alternans).
# - However, one may in fact leverage this "problem" for high-throughput
# detection of discordant (versus concordant) alternans - whenever the
# variable 'alternansData' above differs markedly from the value of
# 'alternansMap' variable, it indicates that alternans is more likely to be
# discordant
