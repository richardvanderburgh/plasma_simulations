%FIRST BLAND cold- Set initial input values

%Input Variables
L = 4*pi; % Physical length of system in meters
nsp = 1; % Number of particle species
nt = 150; % Number of time steps
dt = .1; % Time step in seconds
epsi = 1; % 1 over epsilon naught (F/m) Epsilon normalized to 1
ng = 256; % Number of spatial grid points - only change the power of 2
iw = 2; %1 for nearest grid point, 2 for linear
a1 = 0; %Smoothing factor
a2 = 0;
e0 = 0; % Add uniform electric field e0*cos(w0*time)
w0 = 0;

%Plotting Intervals
irho = 0; %40
irhos = 0;
iphi = 0; %10
ie = 0;
mplot = 3;
ixvx = 1; %5
ivxvy = 1;
ifvx = 0;
nplot = 30; %??


%% Species Input Variables
N = [2048 0]; % Number of simulation particles
wp=[1 1]; %Plasma Frequency
wc = [0 0]; % Cyclotron frequency
qm = [-1 0]; % q/m charge to mass ratio (C/kg)
vt1 = [0 0]; % RMS thermal velocity for random velocities
vt2 = [0 0]; % RMS thermal velocity for ordered velocities
nv2 = [0 0]; 
nlg = [1 0]; %Number of loading groups
v0 = [0 0]; % Drift velocity
pch = [0 0]; % species pitch angle
distribution =0; %Distribution 0=Cold 1=Two-Stream 

%Perturbation
mode=[2 0];
x1=[0.001 0.00];
v1 = [0 0];
thetax=[0 0];
thetav=[0 0];


