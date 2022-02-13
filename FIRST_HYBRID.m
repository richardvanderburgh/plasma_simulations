%FIRST HYBRID- Set initial input values
example='Hybrid Oscillation';
%Input Variables
L = 6.28318530717958; % Physical length of system in meters
nsp = 1; % Number of particle species
nt = 150; % Number of time steps
dt = .1; % Time step in seconds
epsi = 1; % 1 over epsilon naught (F/m) Epsilon normalized to 1
ng = 256; % Number of spatial grid points - only change the power of 2
iw = 1; %1 for nearest grid point, 2 for linear
a1 = 0; %Smoothing factor
a2 = 1000;
e0 = 0; % Add uniform electric field e0*cos(w0*time)
w0 = 0;

%Plotting Intervals
irho = nt/6+1; %40
iphi = nt/6+1; %10
iE = nt/6+1;
mplot = 3;
ixvx = nt/6+1; %60
ivxvy = 0;
ifvx = nt/6+1;
nplot = 30; %??

N = 128; % Number of simulation particles
wp= 1; %Plasma Frequency
wc = -1; % Cyclotron frequency
qm = -1; % q/m charge to mass ratio (C/kg)
vt1 = 0; % RMS thermal velocity for random velocities
vt2 = 0; % RMS thermal velocity for ordered velocities
nv2 = 0;
nlg = 1; %Number of loading groups
v0 = 0; % Drift velocity
distribution =1; %Distribution 0=Cold 1=Two-Stream 

%Perturbation
mode=1;
x1=0;
v1 = 0.02;
thetax=0;
thetav=0;
