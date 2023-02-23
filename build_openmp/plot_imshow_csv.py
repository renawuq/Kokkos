#!/usr/bin/env python

# code adapted from http://matplotlib.org/examples/animation/dynamic_image.html

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def read_csv(filename):
    data = np.loadtxt(filename, delimiter=',')
    return data

fig = plt.figure()

index = 0
prefix = "Kokkos2"
#prefix = "Serial"
filename = "out/"+prefix+"_{0:03d}.csv".format(index)
data = read_csv(filename)
im = plt.imshow(data, cmap=plt.get_cmap('seismic'), animated=True,vmin=-0.2,vmax=0.2)

def updatefig(i):

    filename = "out/"+prefix+"_{0:03d}.csv".format(i)
    data = read_csv(filename)

    im.set_array(data)
    return im,

def inifig():

    filename = "out/"+prefix+"_{0:03d}.csv".format(0)
    data = read_csv(filename)

    im.set_array(data)
    return im,

ani = animation.FuncAnimation(fig, updatefig, np.arange(1, 399), init_func=inifig,interval=100, blit=True)

plt.show()
