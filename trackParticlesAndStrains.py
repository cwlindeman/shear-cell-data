#######################################################################################
# Made by Chloe W L
# Jan 28 2023, cleaned up July 20, 2024
#
# Takes csv files from "analyze particles" in imagej and links them
#
#######################################################################################

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv
import random
import time
import sys
import os
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

plt.rcParams["figure.figsize"] = (5,3)
matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)

def distFunc(x0, y0, x, y):
    delta = np.power(np.power(x - x0, 2.0) + np.power(y - y0, 2.0), 0.5)
    return delta

def findStrainsFromLV(LV, origArea, orig=[[1, 0], [0, 1]]):
    # normalize LV and rotate so it's symmetric
    compression = np.linalg.det(LV)
    fractionalCompression = np.linalg.det(LV)/origArea
    #LV = np.matmul(LV,np.linalg.inv(orig))/np.sqrt(fractionalCompression)
    LV = LV/np.sqrt(compression)
    rotAng = np.arctan((LV[1,0] - LV[0,1])/(LV[0,0] + LV[1,1]))
    rotMat = np.array([[np.cos(rotAng), np.sin(rotAng)], \
    [-1*np.sin(rotAng), np.cos(rotAng)]])
    LV = np.matmul(rotMat,LV)

    # calculate the magnitude and direction of strain
    gamma = 2*np.sqrt(1 - (2/(LV[0,0] + LV[1,1]))**2.0)
    sign_simple = 1
    if LV[0,0] < LV[1,1]:
        sign_pure = -1
    else:
        sign_pure = 1

    theta = np.arccos(2.0*LV[0,1]*np.sqrt(1.0-(gamma**2.0)/4.0)/gamma)/2.0

    gammaPure = sign_pure*gamma*np.sin(2*theta)/2.0
    gammaSimple = sign_simple*gamma*np.cos(2*theta)/2.0

    #return [gamma, 180.0*theta/np.pi, compression/origArea - 1.0]
    return [gammaSimple, gammaPure, fractionalCompression-1]

def getColorFromOrientation(dx, dy):
    cmap = matplotlib.cm.get_cmap("jet")
    dr = np.array([dx,dy])
    direction = dr/np.linalg.norm(dr)
    rel_vec = np.array([1,1])/np.sqrt(2)
    # relative to 45 degrees
    c = np.abs(np.dot(direction, rel_vec))
    rgba = cmap(c)
    return [rgba[0], rgba[1], rgba[2]]

# the main deal
def findTrajectories(foldername, cycleEndList=[], buckleFlag=False, particles=True, fpc=4):
# if particles = False, will return strains for all frames (and no particle information)
# if particles = True, will only return stroboscopoic information

    tol = 150
    headerLen = 1

    if "Results_green.csv" in os.listdir(foldername) and particles==True:

        # READ IN PARTICLE DATA
        data = [] # first index is slice; next is particle number; next is x or y
        positionList = []
        slice = 1

        with open(foldername + "/Results_green.csv") as csvfile:
            fileData = csv.reader(csvfile, delimiter=',')
            a = 0
            for row in fileData:
                a += 1
                if a == headerLen:
                    for ii in range(len(row)):
                        if str(row[ii]) == "Slice":
                            index_s = np.copy(ii)
                        elif str(row[ii]) == "Area":
                            index_a = np.copy(ii)
                        elif str(row[ii]) == "XM":
                            index_x = np.copy(ii)
                        elif str(row[ii]) == "YM":
                            index_y = np.copy(ii)
                if a > headerLen:
                    if int(row[index_s]) == slice:
                        positionList.append([float(row[index_x]), float(row[index_y]), int(row[index_a])])
                    if int(row[index_s]) != slice:
                        data.append(positionList)
                        positionList = []
                        positionList.append([float(row[index_x]), float(row[index_y]), int(row[index_a])])
                        slice = int(row[index_s])
        data.append(positionList)

        fpc = int((len(data)-1)/5)
        if particles==False:
            fpc = 1
        data = data[::fpc]
        data = np.array(data)
        numFrames = len(data[:,0,0])
        N = len(data[0,:,0])
        print(str(N) + " particles")
        if numFrames != 6:
            print(str(numFrames) + " frames")

        # reorder so that individual particles stay at the same index over
        # all frames
        for ii in range(numFrames-1):
            x = np.array(data[ii+1,:,0])
            y = np.array(data[ii+1,:,1])
            area = np.array(data[ii+1,:,2])
            for jj in range(N):
                x0 = data[ii,jj,0]
                y0 = data[ii,jj,1]
                a0 = data[ii,jj,2]
                a = np.argmin(distFunc(x0*np.ones((len(x))), y0*np.ones((len(y))), x, y))
                minDist = np.min(distFunc(x0*np.ones((len(x))), y0*np.ones((len(y))), x, y))
                if minDist < tol:
                    data[ii+1,jj,:] = [x[a], y[a], area[a]]
                else:
                    print("frame " + str(ii) + ", particle " + str(jj))
                    data[ii+1,jj,:] =  [-1,-1,-1]
                    #indexList.append(-1)


    # READ IN LV DATA
    frames = []
    allx = []
    ally = []

    # read in side boundary positions and fit to lines
    headerLen = 1
    fileName = foldername + '/Results_blue.csv'
    mList1 = []
    bList1 = []
    mList2 = []
    bList2 = []
    with open(fileName,'r') as csvfile:
        fileData = csv.reader(csvfile, delimiter=',')
        a = 0
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        frameNum = 1
        for row in fileData:
            a += 1
            if a == headerLen:
                for ii in range(len(row)):
                    if str(row[ii]) == "Slice":
                        index_s = np.copy(ii)
                    elif str(row[ii]) == "XM":
                        index_x = np.copy(ii)
                    elif str(row[ii]) == "YM":
                        index_y = np.copy(ii)
            elif a > headerLen:
                if int(row[index_s]) != frameNum:
                    frameNum = int(float(row[index_s]))
                    m1, b1 = np.polyfit(x1, y1, 1)
                    mList1.append(m1)
                    bList1.append(b1)
                    m2, b2 = np.polyfit(x2, y2, 1)
                    mList2.append(m2)
                    bList2.append(b2)
                    x1 = []
                    y1 = []
                    x2 = []
                    y2 = []
                allx.append(float(row[index_x]))
                ally.append(float(row[index_y]))
                if float(row[index_x]) < 1000:
                    x1.append(float(row[index_y]))
                    y1.append(float(row[index_x]))
                else:
                    x2.append(float(row[index_y]))
                    y2.append(float(row[index_x]))

    m1, b1 = np.polyfit(x1, y1, 1)
    mList1.append(m1)
    bList1.append(b1)
    m2, b2 = np.polyfit(x2, y2, 1)
    mList2.append(m2)
    bList2.append(b2)

    #numFrames = len(mList1)
    #print(numFrames)

    # read in top/bottom boundary positions and fit to lines
    fileName = foldername + '/Results_red.csv'
    mList3 = []
    bList3 = []
    mList4 = []
    bList4 = []
    with open(fileName,'r') as csvfile:
        fileData = csv.reader(csvfile, delimiter=',')
        a = 0
        x3 = []
        y3 = []
        x4 = []
        y4 = []
        frameNum = 1
        for row in fileData:
            a += 1
            if a == headerLen:
                for ii in range(len(row)):
                    if str(row[ii]) == "Slice":
                        index_s = np.copy(ii)
                    elif str(row[ii]) == "XM":
                        index_x = np.copy(ii)
                    elif str(row[ii]) == "YM":
                        index_y = np.copy(ii)
            if a > headerLen:
                if int(row[index_s]) != frameNum:
                    frameNum = int(float(row[index_s]))
                    m3, b3 = np.polyfit(x3, y3, 1)
                    mList3.append(m3)
                    bList3.append(b3)
                    m4, b4 = np.polyfit(x4, y4, 1)
                    mList4.append(m4)
                    bList4.append(b4)
                    x3 = []
                    y3 = []
                    x4 = []
                    y4 = []
                allx.append(float(row[index_x]))
                ally.append(float(row[index_y]))
                if float(row[index_y]) < 1000:
                    x3.append(float(row[index_x]))
                    y3.append(float(row[index_y]))
                else:
                    x4.append(float(row[index_x]))
                    y4.append(float(row[index_y]))
    m3, b3 = np.polyfit(x3, y3, 1)
    mList3.append(m3)
    bList3.append(b3)
    m4, b4 = np.polyfit(x4, y4, 1)
    mList4.append(m4)
    bList4.append(b4)

    mList1 = np.array(mList1)
    mList2 = np.array(mList2)
    mList3 = np.array(mList3)
    mList4 = np.array(mList4)
    bList1 = np.array(bList1)
    bList2 = np.array(bList2)
    bList3 = np.array(bList3)
    bList4 = np.array(bList4)

    #plt.show()
    #plt.clf()

    # find the center of the box and LV in each frame
    centerx = []
    centery = []
    LVs = []
    for ii in range(len(mList1)):
        # get the slopes at this frame
        m1 = 1/mList1[ii]
        b1 = -1*bList1[ii]/mList1[ii]
        m2 = 1/mList2[ii]
        b2 = -1*bList2[ii]/mList2[ii]
        m3 = mList3[ii]
        b3 = bList3[ii]
        m4 = mList4[ii]
        b4 = bList4[ii]

        # find where all the lines intersect
        i1x = (b3 - b1)/(m1 - m3)
        i2x = (b4 - b1)/(m1 - m4)
        i3x = (b4 - b2)/(m2 - m4)
        i4x = (b3 - b2)/(m2 - m3)
        i1y = m3*i1x + b3
        i2y = m4*i2x + b4
        i3y = m4*i3x + b4
        i4y = m3*i4x + b3

        centerx.append((i1x + i3x)/2.0)
        centery.append((i1y + i3y)/2.0)

        # get the lattice vectors
        LV1a = np.array([i3x - i2x, i3y - i2y])
        LV1b = np.array([i4x - i1x, i4y - i1y])
        LV2a = np.array([i2x - i1x, i2y - i1y])
        LV2b = np.array([i3x - i4x, i3y - i4y])
        LV1 = (LV1a + LV1b)/2.0
        LV2 = (LV2a + LV2b)/2.0
        LVs.append(np.array([[LV1[0], LV2[0]], [LV1[1], LV2[1]]]))


    gamma_s = []
    gamma_p = []
    gamma_c = []
    area = np.linalg.det(LVs[0])
    fpc = int((len(LVs)-1)/5)
    for ii in range(len(LVs)):
        #print(LVs[ii])
        #print(findStrainsFromLV(LVs[ii], area))
        gamma_s.append(findStrainsFromLV(LVs[ii], area, LVs[0])[0])
        gamma_p.append(findStrainsFromLV(LVs[ii], area, LVs[0])[1])
        gamma_c.append(findStrainsFromLV(LVs[ii], area, LVs[0])[2])

    gamma_s = np.array(gamma_s)
    gamma_p = np.array(gamma_p)
    gamma_c = np.array(gamma_c)


    if particles==False:
        fpc = 1
    gamma_s = gamma_s[::fpc]
    gamma_p = gamma_p[::fpc]
    gamma_c = gamma_c[::fpc]
    LVs = LVs[::fpc]


    gammas = np.array([gamma_s, gamma_p, gamma_c])

    if "Results_green.csv" in os.listdir(foldername) and particles==True:

        # subtract off affine motion
        for ii in range(numFrames):
            for jj in range(N):
                data[ii,jj,:2] =  np.dot(np.linalg.inv(LVs[ii]),data[ii,jj,:2])

        # subtract off average motion (drift)
        for ii in range(numFrames):
            for jj in range(N):
                data[ii,:,0] -= np.mean(data[ii,:,0])
                data[ii,:,1] -= np.mean(data[ii,:,1])

        if "readout" in foldername:

            RMS = []
            fpc = 1

            for ii in range(int((numFrames-1)/fpc)+1):
                deltas = []
                for jj in range(N):
                    dd = distFunc(data[fpc*ii,jj,0], data[fpc*ii,jj,1], \
                                      data[0,jj,0], data[0,jj,1])
                    if dd < 0.5:
                        deltas.append(dd)
                RMS.append(np.mean(deltas))

        else:

            RMS = []
            fpc = 1

            for ii in range(int((numFrames - 1)/fpc)):
                deltas = []
                for jj in range(N):
                    deltas.append(distFunc(data[fpc*ii,jj,0], data[fpc*ii,jj,1], \
                                      data[fpc*ii+fpc,jj,0], data[fpc*ii+fpc,jj,1]))
                RMS.append(np.mean(deltas))

    else:
        RMS = []

    return gammas, RMS


# UNCOMMENT ONE LINE TO LOOK AT CORRESPONDING DATA
#fileNames = ["24-06-02-a", "24-06-02-b", "24-06-02-c", "24-06-03-f", "24-06-03-g"] # no training, SB readout
#fileNames = ["24-05-28-a", "24-05-28-b", "24-05-28-c", "24-05-28-d", "24-06-03-h"] # SB training (small), SB readout
#fileNames = ["24-05-28-e", "24-05-28-f", "24-05-28-g", "24-05-28-h", "24-06-03-i"] # SB training (large), SB readout
#fileNames = ["24-06-06-a", "24-06-06-b", "24-06-07-a", "24-06-07-b", "24-06-07-e"] # SA training (small), SA readout
#fileNames = ["24-06-06-c", "24-06-06-d", "24-06-07-c", "24-06-07-d", "24-06-07-f"] # SA training (large), SA readout
#fileNames = ["24-05-29-d", "24-05-29-e", "24-05-29-f", "24-06-03-j", "24-06-03-k"] # SA training, SB readout
fileNames = ["25-03-10-a", "25-03-10-b", "25-03-10-c", "25-03-10-d", "25-03-10-e", "25-03-11-a", "25-03-11-b"] # alternating training, mixed readout

# Will run into trouble with these last set because some experiments were done with stroboscopic-only measurements.
# Can look at deltas only, or else comment out some of the experiments to see gammas.
#fileNames = [ "24-05-29-m", "24-05-29-n", "24-05-29-o", "24-06-02-d", "24-06-03-c"] # mixed training, mixed readout

# GET STRAINS
gammas = []
datas = []
for ii in range(len(fileNames)):
    gamma_r, data_r = findTrajectories(fileNames[ii]+"-readout", particles=False)
    if os.path.isdir(fileNames[ii]):
        gamma_t, data_t = findTrajectories(fileNames[ii], particles=False)
        gammas.append(np.concatenate([np.array(gamma_t),np.array(gamma_r)], axis=1))
    else:
        gammas.append(np.array(gamma_r))

# PLOT STRAINS
gammas = np.mean(gammas,axis=0)
plt.plot(gammas[1,:] - gammas[1,0], '-o', label="S1")
plt.plot(gammas[0,:] - gammas[0,0], '-o', label="S2")
plt.plot(gammas[2,:] - gammas[2,0], '-o', label="C")
plt.legend()
plt.xlabel("frame number")
plt.ylabel(r"$\gamma$")
plt.tight_layout()
plt.show()
plt.clf()

# GET DELTAS
gammas = []
datas = []
for ii in range(len(fileNames)):
    gamma_r, data_r = findTrajectories(fileNames[ii]+"-readout", particles=True)
    datas.append(np.array(data_r))

# PLOT STRAINS
datas = np.mean(datas,axis=0)
plt.plot(datas, '-o')
plt.xlabel("frame number")
plt.ylabel(r"$\Delta_{min}$")
plt.tight_layout()
plt.show()
plt.clf()
