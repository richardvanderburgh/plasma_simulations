% ACCEL - Calculate force and advance velocity

%a = zeros(max(N),nsp);

%p = zeros(max(N),nsp);

if  t==0
    q=-0.5*q;
    ke(t+1,:)=dum;

end

for species=1:nsp
    if t~=0
       q(species) = L*wp(species)*wp(species)/(epsi*N(species)*qm(species));
    end

%  Renormalize acceleration if need be.

      dxdt = dx/dt;
      
      ae = (q(species)/m(species))*(dt/dxdt);
      if T(species)~=0
          ae = 0.5*ae;
      end
      if ae~=ael
       tem = ae/ael;
       for j = 1:ng+1
        a(j) = a(j)*tem;
       end
       ael = ae;
      end
      aedata(t+1,species)=ae;
      aeldata(t+1,species)=ael;
%}
% Redistribute force to particles

% For weight equal to zero (NGP)
if iw == 1   
    %  NGP, grid points at j*dx.
      v1s = 0;
      v2s = 0;
      for i=1:N(species)
       j = floor(x(i) + 0.5)+1;
       vo = vx(i);
       vn = vo + a(j);
       v1s = v1s + vn;
       v2s = v2s + vn*vo;
       v1sdata(t+1,species)=v1s;
       v2sdata(t+1,species)=v2s;
       vx(i) = vn;
      end
      
      v2sdata(t+1,species)=v2s;
      p(t+1,species) = m(species)*v1s*dxdt;
      ke(t+1,species) = 0.5*m(species)*v2s*dxdt^2;
      de(t+1,species)=N(species)*m(species)/2*mean(vx(:,species))^2*dxdt^2; 
      therme(t+1,species)=ke(t+1,species)-de(t+1,species);
%{
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
 %}
    %Linear weighting
elseif iw == 2
    %  Linear, momentum conserving.
     if T(species)==0       
      v1s = 0;
      v2s = 0;
      for i = 1:N(species)
       j = floor(x(i,species));
       ajdata(i,species)=a(j+1);
       jdata(i,species)=j;
       vo = vx(i,species);
       vcont(i,species) = (x(i,species) - floor(x(i,species)))*(a(j+2) - a(j+1));
       vn = vo + a(j+1) + (x(i,species) - j)*(a(j+2) - a(j+1));
       v1s = v1s + vn;
       v1sdata(t+1,species)=v1s;
       v2s = v2s + vo*vn;
       v2sdata(t+1,species)=v2s;
       vx(i,species) = vn;
      end 
    
      p(t+1,species) = m(species)*v1s*dxdt;
      ke(t+1,species) = 0.5*m(species)*v2s*dxdt^2;
      de(t+1,species)=N(species)*m(species)/2*mean(vx(:,species))^2*dxdt^2; 
      therme(t+1,species)=therme(t+1,species)+ke(t+1,species)-de(t+1,species);
     
       
        elseif T(species)~=0 
     %Linear, momentum conserving, uniformly magnetized.
      s = 2.0*T(species)/(1.0 + T(species)^2);
      v2s = 0;
      for i = 1:N(species)
       j = floor(x(i,species))+1;
       aa = a(j) + (x(i,species) - j-1)*(a(j+1) - a(j));
       vyy = vy(i,species);
       vxx = vx(i,species) - t*vyy + aa;
       vyy = vyy + s*vxx;
       vxx = vxx - t*vyy;
       v2s = v2s+ vxx*vxx + vyy*vyy;
       vx(i,species) = vxx + aa;
       vy(i,species) = vyy;
      end
      ke(t+1,species) = ke(t+1,species) + 0.5*m(species)*v2s*dxdt*dxdt;
      de(t+1,species)=N(species)*m(species)/2*mean(vx(:,species))^2*dxdt^2; 
      therme(t+1,species)=therme(t+1,species)+ke(t+1,species)-de(t+1,species);
        end
end
%{
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
        %}
end

%% Kinetic energy and momemtum, drift and thermal energies
%{
for species=1:nsp
    for K = 1:N(species)
        ke(t+1,species)=ke(t+1,species)+m(species)*v_old(K,species)*vx(K,species)/2;
       
        p(K,species)=m(species)*vx(K,species);
    end
    
        de(t+1,species)=N(species)*m(species)/2*mean(vx(:,species))^2;
        
        therme(t+1,species)=therme(t+1,species)+ke(t+1,species)-de(t+1,species);
end       

%}
