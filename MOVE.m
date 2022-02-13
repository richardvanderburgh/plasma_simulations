% MOVE -   Advances position one time step and accumulates charge density.
   
for species=1:nsp
% Clear out old charge density.
for j = 2:ng+1
     rho(j) = rho0(j);
end
rho(1) = 0;
end

for species=1:nsp
    
for i = 1:N(species)   
    
%NGP 
if iw==1
       x(i,species) = x(i,species) + vx(i,species);
       if x(i,species)<0
           x(i,species) = x(i,species) + ng;
       end
       if x(i,species)>ng 
           x(i,species) = x(i,species) - ng;
       end
      j = floor(x(i,species) + 0.5+1);
       rho(j) = rho(j) + qdx(species);
%Linear
elseif iw==2
    
       x(i,species) = x(i,species) + vx(i,species);
       if x(i,species)<0
          x(i,species) = x(i,species) + ng;
       end
       if x(i,species)>=ng 
          x(i,species) = x(i,species) - ng;
       end
    
       j = floor(x(i,species));
       jdata(i,species)=j;
       drho = qdx(species)*(x(i,species) - j);
       drhodata(i,species)=drho;
       rho(j+1) = rho(j+1)-drho+qdx(species);
       rho(j+2) = rho(j+2)+drho;   
       rhojdata(i,species)=rho(j+1);
       
end
 
end
end 


		



%{
for k=1:N(species)
    % determine the bin in which the particle is 
    %fprintf('particle  %d %f %f %d \n', k, dx, x(k), ng);
    jbin = x(k,species)/dx;
    %fprintf(' index %f\n', jbin);
    
    if(jbin < 0); jbin = jbin + ng; end
    if(jbin >= ng); jbin = jbin - ng; end
    
    j = int16(floor(jbin));
    j = j+1; % matlab index start from 1
    %fprintf(' grid piont to add charge %f\n', jbin);
    qj = q(species)*(dx-(x(k,species)-gridx(j)))/dx;
    qjp1 = q(species)*(x(k,species)-gridx(j))/dx;
    rho(j)=rho(j) + qj/dx;
    rho(j+1)=rho(j+1) + qjp1/dx;
end
end
   rho(1) = rho(1)+rho(ng+1);
   rho(ng+1) = rho(1);

   %{
    x_hist(t) = x(1);
    xi_hist(t) = ionx(1);
    v_hist(t) = vx(1);
    vi_hist(t) = vxi(1);
    %}
%}