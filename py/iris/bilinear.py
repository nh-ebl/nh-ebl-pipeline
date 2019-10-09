################################################################################
# Name : mosaiq
# Purpose : This does a transformation on the input image through a
# bilinear interpolation method
#Author : Benjamin Vaughan
#Start Date : Oct 4, 2019
#Additional Info
#
################################################################################
import numpy as np
from scipy.interpolate import griddata as interp
import matplotlib.pyplot as plt

def mbilinear(x, y, array):
    array = np.squeeze(array)
    six = x.shape
    siy = y.shape
    sia = array.shape
    Nax = sia[0]
    Nay = sia[1]
    Nx = six[0]
    Ny = six[1]
    output = np.zeros((Nx, Ny))
    minval = np.min(array)

    arr = np.array([[1, 3, 4],
                    [2, 5, 6]])

    print(arr[:,0])

    count = 0
    for i in range(Nx):
        for j in range(Ny):
            if x[i,j] < 0 or x[i,j] > Nax-1 or y[i,j] < 0 or y[i,j] > Nay-1:
                count += 1

    inter_percent = 1. * (Nx * Ny - count) / (Nx * Ny) * 100
    print('Images Intersection = %s percent' % (inter_percent))

    pixx = np.arange(0, Nax)
    pixy = np.arange(0, Nay)

    for i in range(Ny):
        mind = []
        for j in range(Nx):
            if x[j,i] >= 0 and x[j,i] <= Nax-1 and y[j,i] >= 0 and y[j,i] <= Nay-1:
                mind.append(j)
        mind = np.asarray(mind)
        if len(mind) != 0:
            xx = np.squeeze(x[mind, i])
            yy = np.squeeze(y[mind, i])
            t = interp((pixx, pixy),array, (xx, yy))
            print(t)
            # print(trucx[0])
            output[mind, i] = truc #maybe this needs to be changed.

    #remove values affected by indef values (for highly < 0 indef and generaly > 0 im)
    ind = []
    for i in range(output.shape[0]):
        for j in range(output.shape[1]):
            ind.append([i,j])
    ind = np.asarray(ind)
    output[ind] = np.nan #i guess nan is the same as missing but not sure
    return output
