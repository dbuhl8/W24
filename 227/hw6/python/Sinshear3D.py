import torch
import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from matplotlib import cm


gpu = torch.cuda.is_available()
print(gpu)
device = torch.device("cuda" if gpu else "cpu")


Nmax = 100
kxmin = 0.01
kxmax = 1.5
kzmin = 0.01
kzmax = 1.5

nk = 100
ky = 1.0
Pe = 1000.0
Re = 1000.0
Ri = 0.1

params = "Re" + str(Re) + "Pe" + str(Pe) + "Ri" + str(Ri)

f = 0 # floquet coefficient

#Internal variables for code
maxn = 5*(2*Nmax+1)


kx = torch.tensor(np.zeros(nk))
kz = torch.tensor(np.zeros(nk))
eigenvals = torch.tensor(np.zeros((nk, nk)))

# loop over kx values.
for i in range(0,nk):
    kx[i] = kxmin + (i)*(kxmax-kxmin)/(nk-1)
    for j in range(0, nk):
        kz[j] = kzmin + (j)*(kzmax-kzmin)/(nk-1)
        Amat = torch.tensor(np.zeros((maxn, maxn)))
        Bmat = torch.tensor(np.zeros((maxn, maxn)))
        #realvalues = torch.tensor(np.zeros((maxn, maxn)))
        #D = torch.tensor(np.zeros((maxn, maxn)))
        # enter A and B
        for n in range(-Nmax, Nmax+1):
            indu = n+Nmax
            indv = 2*Nmax+1+indu
            indw = 2*Nmax+1+indv
            indt = 2*Nmax+1+indw
            indp = 2*Nmax+1+indt
            #Dissipative Coefficient
            dispterm = (kx[i]**2 + ky**2*(f+n)**2 + kz[j]**2)
            #Off-diagonal Terms
            if n>-Nmax:
                Amat[indu,indu-1] = -kx[i]/2
                Amat[indu,indv-1] = -ky/2
                Amat[indw,indw-1] = -kx[i]/2
                Amat[indt,indt-1] = -kx[i]/2
                Amat[indv, indv-1] = -kx[i]/2
            if n<Nmax:
                Amat[indu,indu+1] = kx[i]/2
                Amat[indu,indv+1] = -ky/2
                Amat[indw,indw+1] = kx[i]/2
                Amat[indt,indt+1] = kx[i]/2
                Amat[indv,indv+1] = kx[i]/2

            #Dissipitive Terms
            Amat[indu,indu]= -dispterm/Re
            Amat[indw,indw]= -dispterm/Re
            Amat[indv,indv]= -dispterm/Re
            Amat[indt,indt]= -dispterm/Pe

            #Stratification Terms here
            Amat[indt,indw]=-1.0
            Amat[indw,indt]= Ri

            #Pressure Terms
            Amat[indu,indp]=-kx[i]  # avoid complex matrix by defining variable ip
            Amat[indv,indp]=-ky*(f+n)
            Amat[indw,indp]=-kz[j]

            # Eigenvalue Matrix B, 1 on diagonal
            Bmat[indu,indu]=1.0
            Bmat[indv,indv]=1.0
            Bmat[indw,indw]=1.0
            Bmat[indt,indt]=1.0
            Bmat[indp,indp] = 0.0

            # Continuity Equation
            Amat[indp,indu] = kx[i]
            Amat[indp,indw] = kz[j]
            Amat[indp,indv] = ky*(f+n)


        # Solve generalized eigenproblem
        thisval = math.floor(maxn/3) - 1
        Amat = Amat.to(device) 
        Bmat = Bmat.to(device)
        X = torch.tensor(np.zeros((maxn, thisval)), device=device)
        print(X.device)
        V,D = torch.lobpcg(Amat,thisval,Bmat, X)
        # This helps ignore negative and Inf eigenvalues
        realvalues = torch.real(torch.diag(D))
        for k in range(0, thisval):
            if realvalues[k].item()<=0.:
                realvalues[k]= 0.
            if realvalues[k].item()>=1.e2 :
                realvalues[k] = 0
        eigenvals[i, j] = torch.max(realvalues).item();


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

KX, KZ = torch.meshgrid(kx, kz)
surf = ax.plot_surface(KX, KZ, eigenvals[:, :], cmap=cm.coolwarm)
ax.set_xlabel("kx")
ax.set_ylabel("kz")
ax.set_zlabel("Re(Lambda)")
plt.show()

#        f = open("kx.txt", "w")

#        f.close()
        
#        f = open("kx.txt", "w")

#        f.close()

#        f = open("kx.txt", "w")

#        f.close()




