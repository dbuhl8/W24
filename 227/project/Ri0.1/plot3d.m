

close all;
clear;
clc;

Ri = 0.1;

datfileID = fopen('nk.dat', 'r');
nk = fscanf(datfileID, '%d');

kxfileID = fopen('kx.dat','r');
kzfileID = fopen('kz.dat','r');
lfileID = fopen('lambda.dat','r');

kformat = '%f';
ksize = [1 nk];

lformat = '%f';
lsize = [nk nk];


kx = fscanf(kxfileID, kformat, ksize);
kz = fscanf(kzfileID, kformat, ksize);
lambda = fscanf(lfileID, lformat, lsize);


[KX, KZ] = meshgrid(kx, kz);

fig = figure(1);

surf(KX, KZ, lambda);
xlabel("$k_x$", 'Interpreter', 'latex')
ylabel("$k_z$", 'Interpreter', 'latex')
zlabel("Re($lambda$)", 'Interpreter', 'latex')

fig = figure(2);

contourf(KX, KZ, lambda)
xlabel("kx")
ylabel("kz")
title("Contour plot of Re($\lambda$), Ri = " + string(Ri), 'Interpreter', 'latex')
%print("ContourPlotRi="+ string(Ri) + ".pdf", "-dpdf")

fig = figure(3);
hold on;
for i = 1:nk
        plot(kx, lambda(:, i), '-o', 'DisplayName', strcat("kz = ", string(kz(i))))
end
xlabel("kx")
ylabel("Re($\lambda$)", 'Interpreter', 'latex')
legend()
orient('landscape')
print(fig, "shearRi0,1", "-dpdf", '-fillpage')



