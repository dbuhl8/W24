import numpy as np
import matplotlib.pyplot as plt


a = -0.5
b = -4.655
c = 12.809

px = np.array([1, -3, np.pi])
py = np.array([2, 2, np.e])
pz = np.array([3, 5, -np.sqrt(2)])

xstart = -4
xstop = 4

ystart = -4
ystop = 4

step = 0.1

numstep = int((xstop - xstop)/step)

x = np.arange(xstart, xstop, step)
y = np.arange(xstart, xstop, step)

xx, yy = np.meshgrid(x, y)

z = a*xx + b*yy + c

#u = np.linspace(xstart, xstop, 100)
#v = np.linspace(xstart, xstop, 100)
#x = 10 * np.outer(u, v)
#y = 10 * np.outer(, np.sin(v))
#z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))

fig = plt.figure()

ax = fig.add_subplot(projection='3d')
ax.plot_surface(xx, yy, z, linewidth=0)
ax.set_xlim([xstart, xstop])
ax.set_ylim([ystart, ystop])
ax.set_zlim([-10, 10])

ax.scatter(px, py, pz)

plt.show()

plt.savefig('question5.png')
