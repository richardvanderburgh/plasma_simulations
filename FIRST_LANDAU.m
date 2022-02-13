%FIRST LANDAU - Set initial input values
example='Landau Damping ';
%Input Variables
L = 6.28318530717958; % Physical length of system in meters
nsp = 1; % Number of particle species
nt = 150; % Number of time steps
dt = .1; % Time step in seconds
epsi = 1; % 1 over epsilon naught (F/m) Epsilon normalized to 1
ng = 256; % Number of spatial grid points - only change the power of 2
iw = 2; %1 for nearest grid point, 2 for linear
a1 = 0; %Smoothing factor
a2 = 10000;
e0 = 0; % Add uniform electric field e0*cos(w0*time)
w0 = 0;
distribution =1; %Distribution 0=Cold 1=Two-Stream 


%Plotting Intervals
irho = 0;%nt/6+1; %40
iphi = 0;%nt/6+1; %10
iE = 0;%nt/6+1;
mplot = 3;
ixvx = 1;%nt/6+1; %60
ivxvy = 0;
ifvx = 0;%nt/6+1;
nplot = 30; %??


%% Species Input Variables
N = [8192 0]; % Number of simulation particles
wp=[1 0]; %Plasma Frequency
wc = [0 0]; % Cyclotron frequency
qm = [-1 0]; % q/m charge to mass ratio (C/kg)
vt1 = [0.5 0]; % RMS thermal velocity for random velocities
vt2 = [0 0]; % RMS thermal velocity for ordered velocities
nv2 = [0 0];
nlg = [1 0]; %Number of loading groups
v0 = [0 0]; % Drift velocity

%Perturbation
mode=[1 0];
x1=[0.2 0];
v1 = [0 0];
thetax=[0 0];
thetav=[0 0];
