#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 13:48:38 2017

@author: mmag
"""

import pickle
import numpy as np

# NOTE:
# Chose set condition on the flux
val = "0.006Flux"
with open("{}.pickle".format(val), "rb") as f:
    rf = pickle.load(f)

# NOTE:
# Check the density
val = "0.006n"
with open("{}.pickle".format(val), "rb") as f:
    densDict = pickle.load(f)

# NOTE:
# Chose to probe around rho=0.61
key="0.061250504095786654,0.0,0.6999998263277603"
dens = densDict[key]["n"]
var = rf[key]["nRadialFlux"]
std = var.std()
# Calculate the time
time = densDict[key]["time"]
dt = time[1] - time[0]

# NOTE:
# Chose to have condition on 3 standard deviations
condition = 3*std

# Get indices where the values are above the conditions
indicesMeetingCond = np.where(var >= condition)[0]

# Get contiguous indices
contiguousIndices = []
while len(indicesMeetingCond) !=0:
    # -1 in order not to get out of index on last point
    for i in range(len(indicesMeetingCond)-1):
        if indicesMeetingCond[i]+1 != indicesMeetingCond[i+1]:
            break
    contiguous = indicesMeetingCond[0:i+1]
    indicesMeetingCond = indicesMeetingCond[i+1:]
    contiguousIndices.append(contiguous)

# As i ran to len(indicesMeetingCond)-1, we make a special treatment on
# the last index
contiguousIndices[-2] =\
    np.append(contiguousIndices[-2], contiguousIndices[-1][0])
contiguousIndices = contiguousIndices[:-1]
contiguousIndices = tuple(contiguousIndices)

# NOTE:
# Can from contiguousIndices make waiting time pdf and pdf for pulse
# width
contiguousTimes = tuple(indices*dt for indices in contiguousIndices)
waitingTimes = []
pulseWidths  = []
for times in contiguousTimes:
    pulseWidths.append(times[-1]-times[0])
    mid = int(len(times)/2)
    if len(pulseWidths) == 1:
        # Initialize waiting times
        waitingStart = times[mid]
        continue

    waitingEnd = times[mid]
    waitingTimes.append(waitingEnd - waitingStart)
    waitingStart = waitingEnd.copy()

waitingTimes = tuple(waitingTimes)
pulseWidths  = tuple(pulseWidths)
# From contiguousIndices make the width from the max, and pad it with
# 500%
maxLen = 0
for indices in contiguousIndices:
    maxLen = len(indices) if len(indices) > maxLen else maxLen
windowSize = int(maxLen*5.0)

# Transform to slices
slices = []
for indices in contiguousIndices:
    # Find the mid of the indices
    mid = int(len(indices)/2)
    # +1 for symmetry
    curSlice = slice(indices[mid]-windowSize, indices[mid]+windowSize+1)
    # Guard for the beginning
    if curSlice.start >= 0:
        # Guard for the end
        if curSlice.stop <= len(dens):
            slices.append(curSlice)
slices = tuple(slices)

# +1 for symmetry
timeSlice = np.array(range(-windowSize,windowSize+1))*dt

# Here we will capture the density slices
# NOTE: Will have a positive flux if negative density fluctuation is
#       transported inwards (holes)
posDensSlices = []
negDensSlices = []
mid = int((slices[0].stop - slices[0].start)/2)
for curSlice in slices:
    curDens = dens[curSlice]
    if curDens[mid] >= 0:
        posDensSlices.append(dens[curSlice])
    else:
        negDensSlices.append(dens[curSlice])
posDensSlices = tuple(posDensSlices)
negDensSlices = tuple(negDensSlices)

# Calculate the averages
avgPosDens = np.array(posDensSlices).mean(axis=0)
avgNegDens = np.array(negDensSlices).mean(axis=0)

# Save the time traces
from matplotlib.pylab import plt
for nr, curDens in enumerate(posDensSlices):
    fig, ax = plt.subplots()
    ax.plot(timeSlice, curDens)
    ax.grid()
    plt.savefig("posDens{}.png".format(nr))

for nr, curDens in enumerate(negDensSlices):
    fig, ax = plt.subplots()
    ax.plot(timeSlice, curDens)
    ax.grid()
    plt.savefig("negDens{}.png".format(nr))

fig, ax = plt.subplots()
ax.plot(timeSlice, avgPosDens)
ax.grid()
plt.savefig("avgPosDens.png")

fig, ax = plt.subplots()
ax.plot(timeSlice, avgNegDens)
ax.grid()
plt.savefig("avgNegDens.png")
