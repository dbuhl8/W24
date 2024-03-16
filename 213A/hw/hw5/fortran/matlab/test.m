clear;
close all;

tol = 10^(-12);
A = [3.0, 1.0, 0.0; 1.0, 2.0, 1.0; 0.0, 1.0, 1.0]

normdiff = 1000;
norm = 0;

while (normdiff > tol) 

    [Q, R] = qr(A);

    A = R*Q

    normdiff = sqrt(sum(diag(A).^2,"all") - norm)
    norm = sum(diag(A).^2, "all");

end

Af = A;

