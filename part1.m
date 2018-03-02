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
    C.g = 9.80665;                      % metres (32.1740 ft) per s�
    
mn=0.26*C.m_0; %electron mass
Temp = 300; %Given in kelvin
rTime=10000; %run time in timesteps
MTBC = 0.2e-12;
Vleft = 1;%voltage of left side

%thermal velocity
Vth = sqrt(2*C.kb*Temp/mn);

%establish inital electron positions
%working area 200nm x 100nm
workX=200*10^-9;
workY=100*10^-9;

size=1000;
displaySize=10;

X= rand(2,size);
Y= rand(2,size);

%positions initialize
Xpos(1,:)= X(1,:)*workX;
Ypos(1,:)= Y(1,:)*workY;

colour = rand(1,displaySize);

% for normal distribution of velocity
sigma =sqrt(C.kb*Temp/mn);
mu = Vth/sqrt(2);
MBdist = makedist('Normal',mu,sigma);
Xvel = zeros(1,size);
Yvel = zeros(1,size);
for i=1:1:size
    vel = random(MBdist);
    Xvel(1,i) = vel;
    vel = random(MBdist);
    Yvel(1,i) = vel;
end

%set timestep of function
spacStep = 0.01*workY;
dt = spacStep/Vth;
steps = 1000;

%variable change
Xvel(1,:) = Xvel(1,:)*dt;
Yvel(1,:) = Yvel(1,:)*dt;

%percent scatter
Pscat=1-exp(-(dt/MTBC));

MFPcount = zeros(1,size);

Efield = Vleft/workX;
force = Efield*C.q_0;
acceleration = force/mn;
accelVelocity = acceleration*(dt^2);

figure(1)
%main function
for i = 1:1:steps
    
    %accelerate velocities
    Xvel(:,:) = Xvel(:,:) + accelVelocity;
    
    %scattering  
    scattered=rand(1,size);
    scatterCheck = scattered<=Pscat;
    velocity = (Vth/sqrt(2)).*randn(1,size);
    Xvel(scatterCheck) = velocity(scatterCheck)*dt;
    velocity = (Vth/sqrt(2)).*randn(1,size);
    Yvel(scatterCheck) = velocity(scatterCheck)*dt;
    tvelocity = sqrt((Xvel/dt).^2 +(Yvel/dt).^2);
    MFPcount(~scatterCheck) = MFPcount(~scatterCheck)+spacStep;
    
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
    
    %MFP calculation
    MFP = sum(MFPcount)/size;
    
    %meanVel = sqrt(Ysum + Xsum);
    %MTBCavg = MFP/meanVel;
    
    %plotting here
    prevX(i,:) =Xpos(1,:);
    prevY(i,:) =Ypos(1,:);
    for j = 1:1:displaySize
        plot(prevX(:,j),prevY(:,j),'color',[colour(1,j) 0 j/displaySize])
        
        xlim([0 workX])
        ylim([0 workY])
        hold on
        drawnow
    end
    
end