% SETRHO - c  Converts position to computer normalization and accumulates charge density.

qdx=q/dx;
rho = zeros(ng+1,1); % charge density at the spatial grid points
rhos = zeros(ng+1,nsp);
rho0 = zeros(ng+1,nsp);
qjdata=zeros(max(N),nsp);
qjp1data=zeros(max(N),nsp);
drhodata=zeros(max(N),nsp);



for species=1:nsp
% If is first group of particles, then clear out rho.
if species==1  
      
       for j = 2:ng+1
       rho(j) = rho0(j,1);
       end
      rho(1) = 0;   
       xn = ng;
       dxi = 1.0/dx;
end
%Add on fixed neutralizing charge density
%  (not needed when all species are mobile -but harmless.)

    %rhos(:,species)=N(species)*q(species)/L;

     rhos(species)=0;
    rho0(:,species) = rho0(:,species) - rhos(:,species);

    
    for j = 2:ng+1
      rho(j) = rho(j) - rhos(j,species);
    end
     %rho = rho0;
%{
  n2 = 0;
  for k=1:nsp
    n1 = n2;
    n2 = n1 + N(k);
    for m = (n1+1):n2
      i  = floor(x(m)+ 2.0);
      i1 = i+1;
      s2 = (x(m)+ 2.0 - i)*q(k);
      s1 = q(k) - s2;
      rho(i ) = rho(i ) + s1;
      rho(i1) = rho(i1) + s2;
    end
  end
	
  rho(2)    = rho(2) +rho(nxp2) -rho0(2);
  rho(1)    = rho(nxp1);
  rho(nxp2) = rho(2);
  %}
  %NGP
 if iw==1
    for i =1:N(species)
          x(i,species) = x(i,species)*dxi;
       if x(i,species)<0
          x(i,species) = x(i,species) + xn;
       end
       if x(i,species)>xn 
          x(i,species) = x(i,species) - xn;
       end
          j = floor((x(i,species))+1 + 0.5);
          rho(j) = rho(j) + qdx(species);
    end
%Linear
 elseif iw==2
 for i = 1:N(species)
         x(i,species) = x(i,species)*dxi;
       if x(i,species)<0
            x(i,species) = x(i,species) + xn;
       end
       if x(i,species)>xn 
            x(i,species) = x(i,species) - xn;
       end
       j = floor(x(i,species));
      jdata(i,species)=j;
       
       drho = qdx(species)*(x(i,species) - j);
       drhodata(i,species)=drho;
       rho(j+1) = rho(j+1)-drho+qdx(species);
       rho(j+2) = rho(j+2)+drho;   
 end
 end
 
end
 