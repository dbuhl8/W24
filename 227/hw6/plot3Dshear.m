close all
clear all

Ri = 0.5

params = "Re1000Pe1000Ri" + string(Ri)

kx = readmatrix(strcat(params, "kx.dat"));
kz = readmatrix(strcat(params, "kz.dat"));
lambda = readmatrix(strcat(params, "lambda.dat"));

fig = figure(1);

hold on 
%plot a contour or a plot with different kz overlayed
for i=1:numel(kz)
    plot(kx, lambda(:, i), '-o', 'DisplayName', strcat("kz = ", string(kz(i))))
end


xlabel("kx")
ylabel("Re(lambda)")
legend()
orient('landscape')
print(strcat(params, "3Dshear"), "-dpdf", '-fillpage')


fig = figure(2);

[KX, KZ] = meshgrid(kx, kz);
contourf(KX, KZ, lambda)
xlabel("kx")
ylabel("kz")
title("Contour plot of Re($\lambda$), Ri = " + string(Ri), 'Interpreter', 'latex')
print("ContourPlotRi="+ string(Ri) + ".pdf", "-dpdf")


fig = figure(3);

surf(KX, KZ, lambda)
xlabel("kx")
ylabel("kz")
zlabel("Re($lambda$)", 'Interpreter', 'latex')

