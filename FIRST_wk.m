%FIRST COLD- Set initial input values
example='Cold Plasma ';
%Input Variables
L = 256; % Physical length of system in meters
nsp = 1; % Number of particle species
nt = 1024; % Number of time steps
dt = .1; % Time step in seconds
epsi = 1; % 1 over epsilon naught (F/m) Epsilon normalized to 1
ng = 256; % Number of spatial grid points - only change the power of 2
iw = 2; %1 for nearest grid point, 2 for linear
a1 = 0; %Smoothing factor
a2 = 0;
e0 = 0; % Add uniform electric field e0*cos(w0*time)
w0 = 0;
cv = 8; 

%Plotting Intervals
irho = 0;%nt/6+1;
irhos = 0;
iphi =0;%nt/6+1;
ie = 0;
mplot = 3;
ixvx = 0;
ivxvy = 0;
ifvx =0;
iE = 0;%nt/6+1;
%% Species Input Variables

N = 1024; % Number of simulation particles
wp=2; %Plasma Frequency
wc = -1; % Cyclotron frequency
qm = -1; % q/m charge to mass ratio (C/kg)
vt1 = 0; % RMS thermal velocity for random velocities
vt2 = 0; % RMS thermal velocity for ordered velocities
nv2 = 0;
nlg = 1; %Number of loading groups
v0 = 0; % Drift velocity
pch = [0 0]; % species pitch angle
angle = 0; %Magnetic field angle
distribution =1; %Distribution 0=Cold 1=Two-Stream 


%Perturbation
mode=[0 0];
x1=[0 0];
v1 = [0.00 0];
thetax=[0 0];
thetav=[0 0];