clear;
close all;

tol = 10^(-12);
A = [3.0, 1.0, 0.0; 1.0, 2.0, 1.0; 0.0, 1.0, 1.0]

B = [2.0, 1.0, 3.0, 4.0; 1.0, -3.0, 1.0, 5.0; 3.0, 1.0, 6.0, -2.0; 4.0, 5.0, -2.0, -1.0]

normdiff = 1000;
norm = 0;
I = eye(3);
i = 0;
while (normdiff > tol) 

    mu = A(3, 3)

    [Q, R] = qr(A-mu*I);

    A = R*Q + mu*I

    normdiff = sqrt(sum(diag(A).^2,"all") - norm)
    norm = sum(diag(A).^2, "all");
    i = i + 1
end
Af = A

[V, D] = eig(B)



