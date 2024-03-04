% Enter input variables here
Nmax = 20;    %Number of Fourier modes to keep
kxmin = 0.01; %Min value of kx
kxmax = 2;   %Max value of kx
nk = 20;  % number of kx increments
kz = 1; % Assumes Lz = 2pi

f = 0;   % Floquet coefficients 


%Internal variables for code
maxn = 3*(2*Nmax+1);

% loop over kx values.
for i=1:nk
    kx(i) = kxmin + (i-1)*(kxmax-kxmin)/(nk-1);
    Amat = zeros(maxn,maxn);  % zero A matrix
    Bmat = zeros(maxn,maxn);  % zero B matrix
    realvalues = zeros(maxn,maxn);
    D = zeros(maxn,maxn);
    % Enter A and B
    for n = -Nmax:Nmax
        indu = n+Nmax+1;
        indw = 2*Nmax+1+indu;
        indp = 2*Nmax+1+indw;
        term = (kx(i)^2+(f+n)^2);
        % Enter u equations; these have index indu
        if(n>-Nmax) 
            Amat(indu,indu-1) = -kx(i)/2;
            Amat(indu,indw-1) = -kz/2;
        end
        if(n<Nmax) 
            Amat(indu,indu+1) = kx(i)/2;
            Amat(indu,indw+1) = -kz/2;
        end
        Amat(indu,indp)=-kx(i);  % avoid complex matrix by defining variable ip
        Bmat(indu,indu)=1.0;

        % Enter w equations with index indw
        if(n>-Nmax) 
            Amat(indw,indw-1) = -kx(i)/2;
        end
        if(n<Nmax) 
            Amat(indw,indw+1) = kx(i)/2;
        end
        Amat(indw,indp)=-kz*(f+n);
        Bmat(indw,indw)=1.0;

        % Enter cont equations
        Amat(indp,indu) = kx(i);
        Amat(indp,indw) = kz*(f+n);
        Bmat(indp,indp) = 0.0; 

    end
    % Solve generalized eigenproblem
    [V,D] = eig(Amat,Bmat,'qz');
    % This helps ignore negative and Inf eigenvalues
    realvalues = real(diag(D));
    for j=1:maxn
        if(realvalues(j)<=0.)
         realvalues(j)= 0.;   
        end
        if(realvalues(j)>=1.e2)
           realvalues(j)= 0;
        end
    end
    lambda(i) = max(realvalues);

end

plot(kx,lambda,'-o')
xlabel('k_x')
ylabel('Re(\lambda)')

hold on
