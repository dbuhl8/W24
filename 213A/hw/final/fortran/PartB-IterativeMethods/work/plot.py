import numpy as np
import matplotlib.pyplot as plt
import math

d = [2, 5, 10, 100, 1000]

for i in d:
    x1 = np.loadtxt("GJ"+str(i) + "error.dat")
    x2 = np.loadtxt("GS"+str(i) + "error.dat")
    y1 = np.loadtxt("GJ"+str(i) + "iteration.dat")
    y2 = np.loadtxt("GS"+str(i) + "iteration.dat")
   
    plt.figure()
    plt.semilogy(x1, y1, label='Gauss-Jacobi')
    plt.semilogy(x2, y2, label='Gauss-Seidel')
    plt.xlabel("Number of Iterations")
    plt.ylabel("Computed Error")
    plt.legend()
    plt.savefig("D="+str(i)+"plot.png")






