% Enter input variables here
Nmax = 40;    %Number of Fourier modes to keep
kxmin = 0.01; %Min value of kx
kxmax = 1.5;   %Max value of kx
nk = 30;  % number of kx increments
kz = 1; % Assumes Lz = 2pi
Pe = 1000;
Re = 1000;
Ri = 0.3;


f = 0;   % Floquet coefficients 
hold on


%Internal variables for code
maxn = 4*(2*Nmax+1);
%if including Temp: 4*()
%if including V, TEMP: 5*()

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
        %indv = 2*Nmax+1+indu;
        indw = 2*Nmax+1+indu;%+indv
        indt = 2*Nmax+1+indw;
        indp = 2*Nmax+1+indt;
        term = (kx(i)^2+(f+n)^2);
        dispterm = (kx(i)^2 + kz^2*(f+n)^2);
        % Enter u equations; these have index indu
        if(n>-Nmax) 
            Amat(indu,indu-1) = -kx(i)/2;
            Amat(indu,indw-1) = -kz/2;
            %temp and v equations
            Amat(indt,indt-1) = -kx(i)/2;
            %Amat(indv, indv-1) = -kx(i)/2;
        end
        if(n<Nmax) 
            Amat(indu,indu+1) = kx(i)/2;
            Amat(indu,indw+1) = -kz/2;
            %temp and v equations
            Amat(indt,indt+1) = kx(i)/2;
            %Amat(indv, indv-1) = -kx(i)/2;
        end
        %dissipative terms on the diagonal
        Amat(indu,indu)= -dispterm/Re;
        Amat(indw,indw)= -dispterm/Re;
        %Amat(indv,indv)= -dispterm/Re;
        Amat(indt,indt)= -dispterm/Pe;
        %Stratification Terms here
        Amat(indt,indw)=-1.0;
        Amat(indw,indt)= Ri;

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
        Bmat(indt,indt)=1.0;

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

plot(kx,lambda,'-o','DisplayName', strcat("Ri = ", string(Ri)))
xlabel('k_x')
ylabel('Re(\lambda)')
legend()

