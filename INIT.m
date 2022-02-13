% INIT - Initial values for each species
cv = 8; 
angle = 0; %Magnetic field angle
gridx  = 0:L/ng:L; % spatial grid (bin width)
dx = gridx(2);  % spatial grid bin width
dtdx=dt/dx;
nxp1=ng+1;
nxp2=ng+2;

%renorm
%cs = cv^2;

  slx = ng;
  npt = sum(N(1:nsp));
  nxp1 = ng+1;
  nxp2 = ng+2;
  X1 = 1:ng;
  X2 = 2:(ng+1);
  X3 = 3:(ng+2);
  cs = cv*cv;
  tcs = 2.0*cs;
  
  q = ng./N(1:nsp).*(wp(1:nsp).^2)./qm(1:nsp);
  mass = q./qm(1:nsp);
 % rho0 = -sum(q(1:nsp).*N(1:nsp))/ng;
 % rho0 = rho0*ones(nxp2,1);

  theta = pi/180*angle;
  costh = cos(theta);
  sinth = sin(theta);

  %
  b0 = wc(1)/qm(1);
  bx0 = b0*costh;
  by0 = b0*sinth;
  
E = zeros(ng+1,nt+1);   

%-- Field Initialization --
  ex = zeros(nxp2,1);
  ey = zeros(nxp2,1);
  ez = zeros(nxp2,1);
  by = ones(nxp2,1)*by0;
  bz = zeros(nxp2,1);

  ajx = zeros(nxp2,1);
  ajy = zeros(nxp2,1);
  ajz = zeros(nxp2,1);

%  rho = zeros(nxp2,1);


x=zeros(max(N),nsp);
vx=zeros(max(N),nsp);
vy=zeros(max(N),nsp);
vz=zeros(max(N),nsp);

%% Energies
% Total time for plot

ESE=zeros(nt+1,1);
esestot=zeros(nt+1,1);
esem=zeros(nt+1,mplot+1);

ke=zeros(nt+1,nsp);
p = zeros(nt+1,nsp);
de=zeros(nt+1,nsp);
therme=zeros(nt+1,nsp);
te=zeros(nt+1,1);

%%
%{
v_old=zeros(max(N),nsp);
vy_old=zeros(max(N),nsp);
%}
ddx=zeros(nsp,1);

q=zeros(1,nsp);
for species=1:nsp    

ngr=N(species)./nlg(species);

q(species) = L*wp(species)*wp(species)/(epsi*N(species)*qm(species));
m(species) = q(species)/qm(species);
nm = N(species)*m(species);
T(species) = tan(-wc(species)*dt/2.);
ddx(species)=L/N(species); % N charges in lenght L, charge cloud width
lg = L/nlg(species);

%%%%%%%%%Computer defined variable


% Set evenly spaced charge distribution distribution
% remember that ddx is the width of the charge cloud

% charge location
%x = transpose(linspace(0,L,N));
%ionx = transpose(linspace(0,L,N));

%for I=1:N(species)
%    x(I,species)=(I)*ddx(species);
%end

for I=1:N(species)
    x(I,species)=(I-0.5)*ddx(species);
end

vx(:,species) = v0(species);

%Load ordered velocities in vx ("quiet start", or at least subdued).
% Is set up for Maxwellian*v**nv2, but can do any smooth distribution.
% Hereafter, ngr is preferably a power of 2.
% First store indefinite integral of distribution function in x array.
% Use midpoint rule -simple and quite accurate.

if vt2(species)~=0
vmax=5*vt2(species);
dv=2*vmax/(N(species)-1);
vvnv2=1;

x(1,species) = 0;

for ith = 2:N(species)
vv=((ith-1.5)*dv-vmax)/vt2(species);

if nv2(species) ~=0
    vvnv2=vv^nv2(species);
end
fv=vvnv2*exp(-0.5*vv^2);
 i1 = ith - 1 + 1;
 x(i1,species) = x(i1-1,species) + max(fv,0);
end

% For evenly spaced (half-integer multiples) values of the integral,
% find corresponding velocity by inverse linear interpolation.

df = x(i1,species)/ngr;
i1 = 1;  
j = 1; 
for ith = 0:ngr-1
       fv = (ith+1 - .5)*df;
             
     while(fv >= x(j+1,species)) 
     j = j + 1;
     end
     
      vv = dv*(j-1 - 1 + (fv - x(j-1,species))/(x(j,species) - x(j-1,species))) - vmax;
    vx(i1,species) = vx(i1,species) + vv;
     i1 = i1 + 1;
end
% For ordered velocities, scramble positions to reduce correlations.
% Scrambling is done by bit-reversed counter -compare sorter in cpft.
% xs=.000,.100,.010,.110,.001,.101,.011,.111,.0001.. (binary fractions).
       xs = 0;
 for  ith = 1:ngr
         i1 = ith - 1 + 1;
        x(i1,species) = xs*lg + .5*ddx(species);
         xsi = 1.;
        while xs>=0
        xsi = .5*xsi;
          xs = xs - xsi;
        end
        xs = xs + 2.*xsi;
  
 i1 = ngr + 1 - 1;
 end

end
  %If magnetized, rotate (vx,0) into (vx,vy).
  if wc(species)~=0 
       for  ith = 1:ngr
         i1 = ith - 1 + 1;
         vv = vx(i1,species);
         theta = 2*pi*x(i1,species)/lg;
          vx(i1,species) = vv*cos(theta);
          vy(i1,species) = vv*sin(theta);
       end
  end
     %Copy first group into rest of groups.
   if nlg(species)~=1  
         j = ngr + 1;
       xs = 0;
       for ith = j:N(species)
        xs = xs + lg;
        for j = 1:ngr
         i1 = j - 1 + 1;
        i2 = i1 + ith - 1;
        x(i2,species) = x(i1,species) + xs;
          vx(i2,species) = vx(i1,species);
         if wc(species)~=0
          vy(i2,species) = vy(i1,species);
         end
        end
       end
   end
   
    
      %From ES1.f
      %Add random maxwellian.
if vt1(species)~=0 
      for I = 1:N(species)
      rm = 0;
      for ith = 1:12
      rm = rm + rand;
      end
       rm = rm - 6;
      
      

         
         
         if wc(species)~=0
        vy(I,species) = vy(I,species) + vt1(species)*rm;
         end
        vx(I,species) = vx(I,species) + vt1(species)*rm;
         
      end
end
    
      %From PIC.py
      %Add random maxwellian.
      %{
     if vt1(species)~=0 
         Fv=zeros(N(species),nsp);
         for i=1:N(species)
            fmax = 0.5 * (1. + exp(-2 * v0(species)^2));
            vmin = -5 * v0(species);
            vmax = 5 * v0(species);
            vtemp = vmin + (vmax - vmin)*(rand);
            f = 0.5 * (exp(-(vtemp - v0(species))^2 / 2.0) + exp(-(vtemp + v0(species))^2 / 2.0));
            xtemporal = fmax * (rand);
            while xtemporal>f
                fmax = 0.5 * (1. + exp(-2. * v0(species)^2));
                vmin = -5 * v0(species);
                vmax = 5 * v0(species);
                vtemp = vmin + (vmax - vmin)*(rand);
            f = 0.5 * (exp(-(vtemp - v0(species))^2 / 2.0) + exp(-(vtemp + v0(species))^2 / 2.0));
                Fv(i,nsp) = f;
                xtemporal = fmax * (rand);
            end
            vx(i,species) = vtemp;
         end
     end
      %}
    
for a = 1:N(species)
    theta = 2*(pi)*mode(species)*x(a,species)/L;
    x(a,species)=x(a,species)+x1(species)*cos(theta + thetax(species));
    
    
    vx(a,species)=vx(a,species)+v1(species)*sin(theta + thetav(species));
    
  
     if (x(a,species) >= L)
            x(a,species) = x(a,species)-L;
     end
     if (x(a,species) < 0)
            x(a,species) = x(a,species)+L;
     end
end


end
N_grid = linspace(0,L,32);

gridt = linspace(0,(nt+1)*dt,(nt+1));
x_hist = zeros(nt,1);

xvxt=zeros(nt,1);
 rhot=zeros(nt,1);
 phit=zeros(nt,1);