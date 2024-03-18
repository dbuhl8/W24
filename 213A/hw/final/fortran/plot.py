import matplotlib.pyplot as plt
import numpy as np
import math

#file10 = open("Image_appn_100010.dat", "r")
#file20 = open("Image_appn_100020.dat", "r")  
#file40 = open("Image_appn_100040.dat", "r")
#file80 = open("Image_appn_100080.dat", "r")
#file160 = open("Image_appn_100160.dat", "r")
#file320 = open("Image_appn_100320.dat", "r")
#file640 = open("Image_appn_100640.dat", "r")
#file1279 = open("Image_appn_101279.dat", "r")

A10 = np.loadtxt("Image_appn_100010.dat")
A20 = np.loadtxt("Image_appn_100020.dat")
A40 = np.loadtxt("Image_appn_100040.dat")
A80 = np.loadtxt("Image_appn_100080.dat")
A160 = np.loadtxt("Image_appn_100160.dat")
A320 = np.loadtxt("Image_appn_100320.dat")
A640 = np.loadtxt("Image_appn_100640.dat")
A1279 = np.loadtxt("Image_appn_101279.dat")

plt.figure()
plt.imshow(A10, cmap='gray', vmin=0, vmax=255)
plt.title("Reduced Rank Reconstruction for K = 10")

plt.figure()
plt.imshow(A20, cmap='gray', vmin=0, vmax=255)
plt.title("Reduced Rank Reconstruction for K = 20")

plt.figure()
plt.imshow(A40, cmap='gray', vmin=0, vmax=255)
plt.title("Reduced Rank Reconstruction for K = 40")

plt.figure()
plt.imshow(A80, cmap='gray', vmin=0, vmax=255)
plt.title("Reduced Rank Reconstruction for K = 80")

plt.figure()
plt.imshow(A160, cmap='gray', vmin=0, vmax=255)
plt.title("Reduced Rank Reconstruction for K = 160")

plt.figure()
plt.imshow(A320, cmap='gray', vmin=0, vmax=255)
plt.title("Reduced Rank Reconstruction for K = 320")

plt.figure()
plt.imshow(A640, cmap='gray', vmin=0, vmax=255)
plt.title("Reduced Rank Reconstruction for K = 640")

plt.figure()
plt.imshow(A1279, cmap='gray', vmin=0, vmax=255)
plt.title("Reduced Rank Reconstruction for K = 1279")

plt.show()




