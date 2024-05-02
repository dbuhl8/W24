clear all
clf reset
Nmax = 10;    %Number of Fourier modes to keep
Mmax = 5;
nk = 5;  % number of kx increments
%kxmin = 0.0; %Min value of kx
%kxmax = 1.5;   %Max value of kx
kx = 0.5;
ky = 1;
kzmin = 0.01; %Min value of kz
kzmax = 1.5; %Max value of kz
maxn = 5*(2*Nmax+1)*(2*Mmax+1);
lambda = zeros(nk,1);

%Hello to Nate from Arstan with love

f = 0.0;
Re = 1000;
Pe = 1000;
Ri = 100;

%+1/0/-1
v_l = 0.5;
v_0 = -0.915;
v_r = -0.5;
%Amplitude of u prime 
A = 0.1;

for j=1:nk
    kz(j) = kzmin+(j-1)*(kzmax-kzmin)/(nk-1);
    Amat = zeros(maxn,maxn);  % zero A matrix
    Bmat = zeros(maxn,maxn);  % zero B matrix
    realvalues = zeros(maxn,maxn);
    D = zeros(maxn,maxn);
    for k = (-Mmax):Mmax
        for i = (-Nmax):Nmax
            indu = i+Nmax+1 + (2*Nmax+1)*(Mmax+k);
            indv = (2*Nmax+1)*(2*Mmax+1)+indu;
            indw = (2*Nmax+1)*(2*Mmax+1)+indv;
            indt = (2*Nmax+1)*(2*Mmax+1)+indw;
            indp = (2*Nmax+1)*(2*Mmax+1)+indt;
            term = (kx^2*(f+i)^2 +ky^2*(f+k)^2+kz(j)^2);
            %u-equation here
            if (i>-Nmax) 
                Amat(indu,indu-1) = complex(0,-0.5*(f+k)*ky*A*v_0);
                Amat(indv,indu-1) = complex(0,-0.5*kx*A*v_0);
                Amat(indv,indv-1) = complex(0,-0.5*ky*(f+k)*A*v_0);
                Amat(indw,indw-1) = complex(0,-0.5*(f+k)*ky*A*v_0);
                Amat(indt, indt-1) = complex(0,-0.5*(f+k)*ky*A*v_0);

              if (k>-Mmax)
                Amat(indu,indu-2*Nmax-2) =complex(0, -(kx*A/2)-(f+k-1)*ky*A*v_r/2 - (f+i-1)*kx*A/2);
                Amat(indu, indv-2*Nmax-2) =complex(0, -ky*A/2);
                Amat(indv, indv-2*Nmax-2) =complex(0, -(i+f-1)*kx*A/2 - ky*A*v_r/2 - ((f+k-1)*ky*A*v_r/2));
                Amat(indv, indu-2*Nmax-2) = complex(0, kx*A*v_r/2);
                Amat(indw, indw-2*Nmax-2) = complex(0,-(i+f-1)*kx*A/2-(f+k-1)*ky*A*v_r/2);
                Amat(indt, indt-2*Nmax-2) = complex(0,-(i+f-1)*kx*A/2-(f+k-1)*ky*A*v_r/2);
              end
              if (k<Mmax)
                Amat(indu, indu+2*Nmax) = complex(0,-kx*A/2+(f+k+1)*ky*A*v_l/2-(f+i-1)*kx*A/2);
                Amat(indu, indv+2*Nmax) = complex(0, ky*A/2);
                Amat(indv, indv+2*Nmax) = complex(0,-(i+f-1)*kx*A/2+(f+k+1)*ky*A*v_l/2 + ky*A*v_l/2);
                Amat(indv, indu+2*Nmax) = complex(0, kx*A*v_l/2);
                Amat(indw, indw+2*Nmax) = complex(0,-(i+f-1)*kx*A/2 + (f+k+1)*ky*A*v_l/2);
                Amat(indt, indt+2*Nmax) = complex(0,-(i+f-1)*kx*A/2 + (f+k+1)*ky*A*v_l/2);
              end
            end


            if (i<Nmax)
                Amat(indu, indu+1) = complex(0,-(f+k)*ky*A*v_0/2);
                Amat(indv, indv+1) = complex(0,-(f+k)*ky*A*v_0/2);
                Amat(indv, indu+1) = complex(0, kx*A*v_0*0.5);
                Amat(indw, indw+1) = complex(0,-(f+k)*ky*A*v_0/2);
                Amat(indt, indt+1) = complex(0,-(f+k)*ky*A*v_0/2);
                if (k>-Mmax)
                    Amat(indu, indu-2*Nmax) = complex(0,kx*A/2 - (f+i+1)*kx*A/2 + (f+k-1)*ky*A*v_r/2);
                    Amat(indu, indv-2*Nmax) = complex(0,-ky*A/2);
                    Amat(indv, indv-2*Nmax) = complex(0,-ky*A*v_l/2 - (f+i+1)*kx*A/2 - (f+k-1)*ky*A*v_l/2);
                    Amat(indv, indu-2*Nmax) = complex(0,-kx*A*v_l/2);
                    Amat(indw, indw-2*Nmax) = complex(0,-(i+f+1)*kx*A/2 + (f+k-1)*ky*A*v_l/2);
                    Amat(indt, indt-2*Nmax) = complex(0,-(i+f+1)*kx*A/2 + (f+k-1)*ky*A*v_l/2);
                end
                if (k<Mmax)
                    Amat(indu, indu+2*Nmax+2) = complex(0,kx*A/2 - (f+i+1)*kx*A/2- (f+k+1)*ky*A*v_r/2);
                    Amat(indu, indv+2*Nmax+2) = complex(0,ky*A/2);
                    Amat(indv, indu+2*Nmax+2) = complex(0,kx*A*v_l/2);
                    Amat(indv, indv+2*Nmax+2) = complex(0,-ky*A*v_l/2 - (f+i+1)*kx*A/2 + (f+k+1)*ky*A*v_l/2);
                    Amat(indw, indw+2*Nmax+2) = complex(0,-(i+f+1)*kx*A/2 - (f+k+1)*ky*A*v_l/2);
                    Amat(indt, indt+2*Nmax+2) = complex(0,-(i+f+1)*kx*A/2 - (f+k+1)*ky*A*v_l/2);
                end
            end

            if (k>-Mmax)
                Amat(indu, indu-(2*Nmax+1)) = complex(-(f+i)*kx/2,0);
                Amat(indu, indv-(2*Nmax+1)) = complex(-ky/2,0);
                Amat(indv, indv-(2*Nmax+1)) = complex(-(f+i)*kx/2,0);
                Amat(indw, indw-(2*Nmax+1)) = complex(-(f+i)*kx/2,0);
                Amat(indt, indt-(2*Nmax+1)) = complex(-(f+i)*kx/2,0);
            end

            if (k<Mmax)
                Amat(indu, indu+(2*Nmax+1)) = complex((f+i)*kx/2,0);
                Amat(indu, indv+(2*Nmax+1)) = complex(-ky/2,0);
                Amat(indv, indv+(2*Nmax+1)) = complex((f+i)*kx/2,0);
                Amat(indw, indw+(2*Nmax+1)) = complex((f+i)*kx/2,0);
                Amat(indt, indt+(2*Nmax+1)) = complex((f+i)*kx/2,0);
            end
             
         

            
            %Dissipation terms
            Amat(indu,indu)= -1/Re*term;
            Amat(indv,indv)= -1/Re*term;
            Amat(indw,indw)= -1/Re*term;
            Amat(indt,indt)= -1/Pe*term;

            %Pressure terms
            Amat(indu,indp)= -kx*(f+i);
            Amat(indv,indp)= -ky*(f+k);
            Amat(indw,indp)= -kz(j);

            %Stratification
            Amat(indt,indw)= -1.0;
            Amat(indw,indt)= Ri;
     

            %ICont equation
            Amat(indp,indu)= kx*(f+i);
            Amat(indp,indv)= ky*(f+k);
            Amat(indp,indw)= kz(j);
            
            %B-mat
            Bmat(indu,indu)= 1.0;
            Bmat(indv,indv)= 1.0;
            Bmat(indw,indw)= 1.0;
            Bmat(indt,indt)= 1.0;
            Bmat(indp,indp)= 0.0;

   
    end

    end
[V,D] = eig(Amat,Bmat,"qz");
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
    lambda(j,1) = max(realvalues);
    

    for i=1:maxn
        if (realvalues(i)>=0.25)
            display(i)
        end
    end

end
