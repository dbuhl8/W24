import numpy as np
import matplotlib.pyplot as plt


b = 2.080
a = 1.920

step = 0.001

numstep = int((b - a)/step)

x = np.linspace(a, b, numstep)

f = (x - 2)**9
g = x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 + 5376*x**3 - 4608*x**2 + 2304*x - 512

fig = plt.figure()
fig.suptitle('Question 13a')
ax = fig.add_subplot()
ax.plot(x, f, label='f(x)')
ax.plot(x, g, label='g(x)')
ax.legend()

plt.savefig('question13.png')
