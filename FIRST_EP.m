%FIRST BEAM-PLASMA - Set initial input values
example='Electron Beam Plasma';

%% General Input Variables
L = 2*pi; % Physical length of system in meters
nsp = 2; % Number of particle species
nt = 300; % Number of time steps
dt = .2; % Time step in seconds
epsi = 1; % 1 over epsilon naught (F/m) Epsilon normalized to 1
ng = 64; % Number of spatial grid points - only change the power of 2
iw = 2; %1 for nearest grid point, 2 for linear
a1 = 0; %Smoothing factor
a2 = 0;
e0 = 0; % Add uniform electric field e0*cos(w0*time)
w0 = 0;

%Plotting Intervals
irho = 100;
irhos = 100;
iphi = 50;
iE = 0;
mplot = 4;
ixvx = 1;
ivxvy =1;
ifvx = 100;
nplot = 30; %??

%% Species Input Variables
N = [600, 64]; % Number of simulation particles
wp= [1, 1]; %Plasma Frequency
wc = [0, 0]; % Cyclotron frequency
qm = [-1, 0.001, ]; % q/m charge to mass ratio (C/kg)
vt1 = [0, 0]; % RMS thermal velocity for random velocities
vt2 = [1 0]; % RMS thermal velocity for ordered velocities
nv2 = [0 0]; 
nlg = [1 1]; %Number of loading groups
v0 = [0 0]; % Drift velocity
distribution =2; %Distribution 0=Cold 1=Two-Stream 2=electron plasma

%Perturbation
mode=[0 0];
x1=[0.0 0];
v1 = [0 0];
thetax=[0 0];
thetav=[0 0];
