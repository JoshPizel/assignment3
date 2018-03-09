%clear all

clearvars
clearvars -GLOBAL
close all

global C
global X Y

    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    
mn=0.26*C.m_0; %electron mass
Temp = 300; %Given in kelvin
rTime=10000; %run time in timesteps
MTBC = 0.2e-12;
Vleft = 0.1;%voltage of left side
electronConc = 10e15;

%thermal velocity
Vth = sqrt(2*C.kb*Temp/mn);

%establish inital electron positions
%working area 200nm x 100nm
workX=200*10^-9;
workY=100*10^-9;
area =workX*workY;

size=1000;
displaySize=10;

X= rand(2,size);
Y= rand(2,size);

%positions initialize
Xpos(1,:)= X(1,:)*workX;
Ypos(1,:)= Y(1,:)*workY;

colour = rand(1,displaySize);

% for normal distribution of velocity
Vthn = Vth/sqrt(2);
Xvel = Vthn*randn(1,size);
Yvel = Vthn*randn(1,size);

%set timestep of function
spacStep = 0.01*workY;
dt = spacStep/Vth;
steps = 1000;

%variable change
Xvel(1,:) = Xvel(1,:)*dt;
Yvel(1,:) = Yvel(1,:)*dt;

%percent scatter
Pscat=1-exp(-(dt/MTBC));

%electric field original
Efield = Vleft/workX;
force = Efield*C.q_0;
acceleration = force/mn;
accelVelocity = acceleration*(dt^2);
    
%set boundary conditions
box1X = [

centreX = workX/2;
centreY = workY/2;

%matrices
G = zeros(L*W,L*W);
B = zeros(L*W,1);

%conductivity
s1 = 1;
s2 = 0.01;

%resistive regions size
rL = L*1/4;
rW = W*2/5;

Smap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n = j+(i-1)*W;
        nxm = j+(i-2)*W;
        nxp = j+i*W;
        nyp = j+1+ (i-1)*W;
        nym = j-1+ (i-1)*W;
        
        if(i==1)
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
            Smap(i,j) = s1;
        elseif(i==L)
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
            Smap(i,j) = s1;
        elseif(j==1)
            G(n,:) = 0;
            G(n,n) = -3;
            if(i > centreX-(rL/2) && i < centreX+(rL/2))%x boundary
                G(n,nxm) = s2;
                G(n,nxp) = s2;
                G(n,nyp) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1;
                G(n,nxp) = s1;
                G(n,nyp) = s1;
                Smap(i,j) = s1;
            end
        elseif(j==W)
            G(n,:) = 0;
            G(n,n) = -3;
            if(i > centreX-(rL/2) && i < centreX+(rL/2))%x boundary
                G(n,nxm) = s2;
                G(n,nxp) = s2;
                G(n,nym) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1;
                G(n,nxp) = s1;
                G(n,nym) = s1;
                Smap(i,j) = s1;
            end
        else          
            G(n,:) = 0;
            G(n,n) = -4;
            %  |           X Boundaries                | 
            if((i > centreX-(rL/2) && i < centreX+(rL/2)) && ...
                    (j > centreY+(rW/2) || j < centreY-(rW/2)))
            %       |            Y Boundaries                |
                G(n,nxp) = s2;
                G(n,nxm) = s2;
                G(n,nyp) = s2;
                G(n,nym) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxp) = s1;
                G(n,nxm) = s1;
                G(n,nyp) = s1;
                G(n,nym) = s1;
                Smap(i,j) = s1;
            end     
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

[Ey,Ex] = gradient(Vmap);

E = gradient(Vmap);

J = Smap.*E;

figure(1)
currentHistory = zeros(1,steps);
%Electron positions plotting and calculating
for i = 1:1:steps
    
    %accelerate velocities
    Xvel(:,:) = Xvel(:,:) + accelVelocity;
    
    %scattering  
    scattered=rand(1,size);
    scatterCheck = scattered<=Pscat;
    velocity = Vthn*randn(1,size);
    Xvel(scatterCheck) = velocity(scatterCheck)*dt;
    velocity = Vthn*randn(1,size);
    Yvel(scatterCheck) = velocity(scatterCheck)*dt;
    tvelocity = sqrt((Xvel/dt).^2 +(Yvel/dt).^2);
    
    %position advance
    %logical indexing
    checkXright = Xpos +Xvel>2e-7;%right side period
    Xpos(checkXright) = Xpos(checkXright)+Xvel(checkXright)-workX;
    checkXleft = Xpos +Xvel<0;%left side period
    Xpos(checkXleft) = Xpos(checkXleft) +Xvel(checkXleft)+workX;
    
    %leftover x 
    leftover = ~(checkXright | checkXleft);
    
    Xpos(leftover) = Xpos(leftover) +Xvel(leftover);
    
    %reflect Y boundary
    checkY = (Ypos+Yvel>1e-7 | Ypos+Yvel<0);
    Yvel(checkY) = Yvel(checkY).*(-1);
    Ypos(1,:) = Ypos(1,:)+Yvel(1,:);
    
    %temperature calculations
    calcTemp = 0.5*mn*(tvelocity.^2)/(2*C.kb);
    averageTemp = sum(calcTemp)/size;
    
    
    %current tracking
    avgVel=sum(tvelocity)/size;
    mu = (avgVel)/Efield;
    currentHistory(i) =C.q_0*electronConc*mu*Efield/area;
    
    %plotting here
    prevX(i,:) =Xpos(1,:);
    prevY(i,:) =Ypos(1,:);
%     for j = 1:1:displaySize
%         plot(prevX(:,j),prevY(:,j),'color',[colour(1,j) 0 j/displaySize])
%         
%         xlim([0 workX])
%         ylim([0 workY])
%         hold on
%         drawnow
%     end
end

for j = 1:1:displaySize
        plot(prevX(:,j),prevY(:,j),'color',[colour(1,j) 0 j/displaySize])
        
        xlim([0 workX])
        ylim([0 workY])
        hold on
        drawnow
end

figure(5)
surf(Vmap)
colorbar
title('Voltage map'),xlabel('X'),ylabel('Y'),zlabel('Voltage')

figure(7)
surf(Ex)
colorbar
title('Electric field X'),xlabel('X'),ylabel('Y'),zlabel('E Field');

figure(8)
surf(Ey)
colorbar
title('Electric field Y'),xlabel('X'),ylabel('Y'),zlabel('E Field');

figure(9)
surf(J)
colorbar
title('Current density'),xlabel('X'),ylabel('Y'),zlabel('Current/m^2');

