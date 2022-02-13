%FIRST EE - Set initial input values
example='Electron-Electron Stream ';
%Input Variables
L = 6.28318530717958; % Physical length of system in meters
nsp = 2; % Number of particle species
nt = 600; % Number of time steps
dt = .1; % Time step in seconds
epsi = 1; % 1 over epsilon naught (F/m) Epsilon normalized to 1
ng = 32; % Number of spatial grid points - only change the power of 2
iw = 2; %1 for nearest grid point, 2 for linear
a1 = 0; %Smoothing factor
a2 = 0;
e0 = 0; % Add uniform electric field e0*cos(w0*time)
w0 = 0;

%Plotting Intervals
irho = 0; %nt/10+4;
iphi = 0; %nt/10+4;
iE = 0; %nt/10+4;
mplot = 3;
ixvx = 1; %60
ivxvy = 0;
ifvx = 0;%nt/10+4;
nplot = 30; %??


%% Species Input Variables
N = [1280 1280]; % Number of simulation particles
wp=[1 1]; %Plasma Frequency
wc = [0 0]; % Cyclotron frequency
qm = [-1 -1]; % q/m charge to mass ratio (C/kg)
vt1 = [0 0]; % RMS thermal velocity for random velocities
vt2 = [0 0]; % RMS thermal velocity for ordered velocities
nv2 = [0 0]; 
nlg = [1 1]; %Number of loading groups
v0 = [1 -1]; % Drift velocity
pch = [0 0]; % species pitch angle
distribution =1; %Distribution 0=Cold 1=Two-Stream 


%Perturbation
mode=[1 1];
x1=[0.001 0.001];
v1 = [0 0];
thetax=[0 0];
thetav=[0 0];


