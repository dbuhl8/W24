
params = "Re1000Pe1000Ri0.5"

kx = readmatrix(strcat(params, "kx.dat"));
kz = readmatrix(strcat(params, "kz.dat"));
lambda = readmatrix(strcat(params, "lambda.dat"));

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
surf(KX, KZ, lambda)
xlabel("kx")
ylabel("kz")
zlabel("Re(lambda)")

