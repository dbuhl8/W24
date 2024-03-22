

close all;
clear;
clc;

    datfileID = fopen('nk.dat', 'r');
%    nk = fscanf(datfileID, '%d');
    nk = 20;

    %kxfileID = fopen('kx.dat','r');
    kzfileID = fopen('kz20.dat','r');
    lfileID = fopen('lambda20.dat','r');
    kz2fileID = fopen('kz1.dat','r');
    l2fileID = fopen('lambda1.dat','r');

    kformat = '%f';
    ksize = [1 nk];

    lformat = '%f';
    lsize = [1 nk];


    %kx = fscanf(kxfileID, kformat, ksize);
    kz = fscanf(kzfileID, kformat, ksize);
    lambda = fscanf(lfileID, lformat, lsize);
    kz2 = fscanf(kz2fileID, kformat, ksize);
    lambda2 = fscanf(l2fileID, lformat, lsize);

    deltaL = abs((lambda(1) - lambda2(1))/lambda(1))


    %[KX, KZ] = meshgrid(kx, kz);

    fig = figure(1);

    %surf(KX, KZ, lambda);
    plot(kz, lambda, "DisplayName", "41 Modes")
    hold on;
    plot(kz2, lambda2, "DisplayName", "3 Modes")
    xlabel("$k_z$", 'Interpreter', 'latex')
    ylabel("Re($\lambda$)", 'Interpreter', 'latex')
    title("Re($\lambda$) versus $k_z$", 'Interpreter', 'latex')
    xlim([-0.1, 1.5])
    ylim([.22, 0.3])
    xline(0, '--r', 'DisplayName', '$k_z = 0$')
    legend('Interpreter', 'latex')
    %zlabel("Re($lambda$)", 'Interpreter', 'latex')

%   fig = figure(2);

%   contourf(KX, KZ, lambda)
%   xlabel("kx")
%   ylabel("kz")
%   title("Contour plot of Re($\lambda$), Ri = " + string(Ri), 'Interpreter', 'latex')
%   print("ContourPlotRi="+ string(Ri) + ".pdf", "-dpdf")


