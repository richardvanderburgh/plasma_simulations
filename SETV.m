% SETV - c  Converts particle velocities at t=0 to computer normalization at t=-dt/2.
%{
for species=1:nsp
% Clear out old charge density.
for j = 2:ng+1
     rho(j) = rho0(j);
end
rho(1) = 0;
end
%}
% Rotate v thru angle +0.5*wc*dt and normalize vy.
%  If t=0, there is no magnetic field. the rotation is omitted and no references are made to vy.
for species=1:nsp
if T(species)~=0
       c = 1.0/sqrt(1.0 + T(species)^2);
       s = c*T(species);
       for ii = 1:N(species)
        vxx(:,species) = vx(ii,species);
        vx(ii,species)= c*vxx(:,species) + s*vy(ii,species);
        vy(ii,species) = -s*vxx(:,species) + c*vy(ii,species);
        vy(ii) = vy(ii)*dtdx;
      end
       
end
%Normalize vx.
for K = 1:N(species)
      vx(K,species) = vx(K,species)*dtdx;
end 
end
% Electric impulse to go back 1/2 time step.
dum = 0;

ACCEL
%{
a = zeros(max(N),nsp);
p = zeros(max(N),nsp);

% Redistribute force to particles

for species=1:nsp
    
% For weight equal to zero (NGP)
   
if iw == 1
for j=1:N(species)  
if ((x(j,species) ~= 0) && (x(j,species) ~= 1))      
if (x(j,species) <= (gridx(1)+0.5*gridx(2)))
            a(j,species) = -qm(species)*E(1);
end
end

end

for K=2:ng-1
    for j=1:N(species)
        if ((x(j,species) ~= 0) && (x(j,species) ~= 1)) 
        if (x(j,species) > (gridx(K)-0.5*gridx(2))) && (x(j) <= (gridx(K)+0.5*gridx(2)))
            a(j,species) = -qm(species)*E(K);
        end
        end
    end
end
    for j=1:N(species)
        if ((x(j,species) ~= 0) && (x(j,species) ~= 1)) 
        if x(j,species) > gridx(ng)-0.5*gridx(2) 
            a(j,species) = -qm(species)*E(ng);
        end
        end
    end
 
 %Linear weighting
elseif iw == 2
    for j= 1:ng
        for K = 1:N(species)
            if ((x(K,species) ~= 0) && (x(K,species) ~= L)) 
            if ((x(K,species) < gridx(j+1)) && (x(K,species) > gridx(j)))
                a(K,species) =  qm(species)*((gridx(j+1)-x(K,species))*E(j)/dx + (x(K,species) - gridx(j))*E(j+1)/dx);
            end
            end
            
        end
    end        
end
    
v_old(:,species) = vx(:,species);
vx(:,species) = v_old(:,species) + a(:,species)*dt;

for K = 1:N(species)
        ke(t+1,species)=ke(t+1,species)+m(species)*v_old(K,species)*vx(K,species)/2;
       
        p(K,species)=m(species)*vx(K,species);
end
        de(t+1,species)=N(species)*m(species)/2*mean(vx(:,species))^2;
        
        therme(t+1,species)=therme(t+1,species)+ke(t+1,species)-de(t+1,species);
end       

%}