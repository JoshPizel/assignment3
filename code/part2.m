
global C

C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;
    

%working area
W = 50;
L = W*3/2;

centreX = L/2;
centreY = W/2;

%matrices
G = zeros(L*W,L*W);
B = zeros(L*W,1);

%conductivity
s1 = 1;
s2 = 0.01;

%resistive regions size
rL = L*1/4;
rW = W*2/5;

%create sigma map
Smap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        if((i > centreX-(rL/2) && i < centreX+(rL/2)) && ...
                    (j > centreY+(rW/2) || j < centreY-(rW/2)))
            Smap(i,j) = s2;
        else
            Smap(i,j) = s1;
        end
    end
end
        


for i =1:1:L
    for j =1:1:W
        n = j+(i-1)*W;
        nxm = j+(i-2)*W;
        nxp = j+i*W;
        nyp = j+1+ (i-1)*W;
        nym = j-1+ (i-1)*W;
        
        if(i==1)
            G(n,:) = 0;
            G(n,n) = Smap(i,j);
            B(n) = 1;
        elseif(i==L)
            G(n,:) = 0;
            G(n,n) = Smap(i,j);
            B(n) = 0;
        elseif(j==1)
            G(n,:) = 0;
            G(n,nxm) = (Smap(i-1,j)+Smap(i,j))/2;
            G(n,nxp) = (Smap(i+1,j)+Smap(i,j))/2;
            G(n,nyp) = (Smap(i,j+1)+Smap(i,j))/2;
            G(n,n) = -(G(n,nxm)+G(n,nxp)+G(n,nyp));
        elseif(j==W)
            G(n,:) = 0;
            G(n,nxm) = (Smap(i-1,j)+Smap(i,j))/2;
            G(n,nxp) = (Smap(i+1,j)+Smap(i,j))/2;
            G(n,nym) = (Smap(i,j-1)+Smap(i,j))/2;
            G(n,n) = -(G(n,nxm)+G(n,nxp)+G(n,nym));
        else          
            G(n,:) = 0;
            G(n,nxm) = (Smap(i-1,j)+Smap(i,j))/2;
            G(n,nxp) = (Smap(i+1,j)+Smap(i,j))/2;
            G(n,nyp) = (Smap(i,j+1)+Smap(i,j))/2;
            G(n,nym) = (Smap(i,j-1)+Smap(i,j))/2;
            G(n,n) = -(G(n,nxm)+G(n,nxp)+G(n,nyp)+G(n,nym));   
        end
    end
end

V = G\B;

%remap
Vmap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n=j+(i-1)*W;
        Vmap(i,j) =V(n);
    end
end

[MX,MY] = meshgrid(1:1:W,1:1:L);
[Ey,Ex] = gradient(Vmap);

figure(5)
surf(Vmap)
colorbar
title('Voltage map'),xlabel('X'),ylabel('Y'),zlabel('Voltage')

figure(6)
quiver(MX,MY,Ex,Ey)
title('Quiver plot of electric field'),xlabel('X'),ylabel('Y')
