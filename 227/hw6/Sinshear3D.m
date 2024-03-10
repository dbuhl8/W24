% Enter input variables here
Nmax = 30;    %Number of Fourier modes to keep
kxmin = 0.01; %Min value of kx
kxmax = 1.5;   %Max value of kx
kzmin = 0.01; %Min value of kx
kzmax = 1.5;   %Max value of kx
nk = 70;  % number of kx/kz increments
ky = 1; % Assumes Lz = 2pi
Pe = 1000;
Re = 1000;
Ri = 0.1;

params = strcat("Re", string(Re), "Pe", string(Pe), "Ri", string(Ri))

f = 0;   % Floquet coefficients 


%Internal variables for code
maxn = 5*(2*Nmax+1);
lambda = zeros(nk, nk);

% loop over kx values.
for i=1:nk
    kx(i) = kxmin + (i-1)*(kxmax-kxmin)/(nk-1);
    for j=1:nk
        kz(j) = kzmin + (j-1)*(kzmax-kzmin)/(nk-1);
        Amat = zeros(maxn,maxn);  % zero A matrix
        Bmat = zeros(maxn,maxn);  % zero B matrix
        realvalues = zeros(maxn,maxn);
        D = zeros(maxn,maxn);
        % Enter A and B
        for n = -Nmax:Nmax
            indu = n+Nmax+1;
            indv = 2*Nmax+1+indu;
            indw = 2*Nmax+1+indv;
            indt = 2*Nmax+1+indw;
            indp = 2*Nmax+1+indt;
            %Dissipative Coefficient
            dispterm = (kx(i)^2 + ky^2*(f+n)^2 + kz(j)^2);

            %Off-diagonal Terms
            if(n>-Nmax) 
                Amat(indu,indu-1) = -kx(i)/2;
                Amat(indu,indv-1) = -ky/2;
                Amat(indw,indw-1) = -kx(i)/2;
                Amat(indt,indt-1) = -kx(i)/2;
                Amat(indv, indv-1) = -kx(i)/2;
            end
            if(n<Nmax) 
                Amat(indu,indu+1) = kx(i)/2;
                Amat(indu,indv+1) = -ky/2;
                Amat(indw,indw+1) = kx(i)/2;
                Amat(indt,indt+1) = kx(i)/2;
                Amat(indv,indv+1) = kx(i)/2;
            end

            %Dissipitive Terms
            Amat(indu,indu)= -dispterm/Re;
            Amat(indw,indw)= -dispterm/Re;
            Amat(indv,indv)= -dispterm/Re;
            Amat(indt,indt)= -dispterm/Pe;

            %Stratification Terms here
            Amat(indt,indw)=-1.0;
            Amat(indw,indt)= Ri;

            %Pressure Terms
            Amat(indu,indp)=-kx(i);  % avoid complex matrix by defining variable ip
            Amat(indv,indp)=-ky*(f+n);
            Amat(indw,indp)=-kz(j);

            % Eigenvalue Matrix B, 1 on diagonal
            Bmat(indu,indu)=1.0;
            Bmat(indv,indv)=1.0;
            Bmat(indw,indw)=1.0;
            Bmat(indt,indt)=1.0;
            Bmat(indp,indp) = 0.0; 

            % Continuity Equation
            Amat(indp,indu) = kx(i);
            Amat(indp,indw) = kz(j);
            Amat(indp,indv) = ky*(f+n);


        end
        % Solve generalized eigenproblem
        [V,D] = eig(Amat,Bmat,'qz');
        % This helps ignore negative and Inf eigenvalues
        realvalues = real(diag(D));
        for k=1:maxn
            if(realvalues(k)<=0.)
            realvalues(k)= 0.;   
            end
            if(realvalues(k)>=1.e2)
            realvalues(k)= 0;
            end
        end
        lambda(i, j) = max(realvalues, [ ], 'all');
    end
end

%need to change this to a 2D plot
[KX, KZ] = meshgrid(kx, kz);
writematrix(kx, strcat(params, "kx.dat"))
writematrix(kz, strcat(params, "kz.dat"))
writematrix(lambda, strcat(params, "lambda.dat"))
%disp(KX)
%disp(KZ)
%disp(lambda)
fig = figure(2);
surf(KX, KZ, lambda)
xlabel('k_x')
ylabel('k_z')
zlabel('Re(lambda)')
%title('Re(lambda)')

hold on
