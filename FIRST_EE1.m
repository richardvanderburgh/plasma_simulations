%FIRST EE - Set initial input values
example='Two Electron Streams';

%Input Variables
L = 6.28318530717958; % Physical length of system in meters
nsp = 2; % Number of particle species
nt = 300; % Number of time steps
dt = .2; % Time step in seconds
epsi = 1; % 1 over epsilon naught (F/m) Epsilon normalized to 1
ng = 32; % Number of spatial grid points - only change the power of 2
iw = 2; %1 for nearest grid point, 2 for linear
a1 = 0; %Smoothing factor
a2 = 0;
e0 = 0; % Add uniform electric field e0*cos(w0*time)
w0 = 0;

%Plotting Intervals
irho = 0;
irhos = 0;
iphi = 0;
iE = 0;
mplot = 3;
ixvx = 1;
ivxvy = 0;
ifvx = 0;

%% Species Input Variables
N = [128 128]; % Number of simulation particles
wp=[1 1]; %Plasma Frequency
wc = [0 0]; % Cyclotron frequency
qm = [-1 -1]; % q/m charge to mass ratio (C/kg)
vt1 = [0 0]; % RMS thermal velocity for random velocities
vt2 = [0 0]; % RMS thermal velocity for ordered velocities
nv2 = [0 0]; 
nlg = [1 1]; %Number of loading groups
v0 = [.2 -.2]; % Drift velocity

%Perturbation
mode=[1 1];
x1=[0.01 -.01];
v1 = [0 0];
thetax=[0 0];
thetav=[0 0];


%% Energies

ESE=zeros(nt+1,1);
esestot=zeros(nt+1);
esesm=zeros(nt+1,mplot);

ke=zeros(nt+1,nsp);
p = zeros(nt+1,nsp);
de=zeros(nt+1,nsp);
therme=zeros(nt+1,nsp);
te=zeros(nt+1,1);

