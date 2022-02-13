% FIELDS - E field solver

rhok=zeros(ng+1,1);
phik=zeros(ng+1,1);
scrach=zeros(ng+1,1);
phi=zeros(ng+1,1);
l=L;
%First time step duties.
%{
      ng2=0;
     if ng2==0  
        ng2 = ng/2;
     end
% Set up ratio phik/rhok.
%  a2>0 gives a short-wavelength cutoff (smoothing).
%  a1>0 gives a mid-range boost to compensate for slight attenuation in force calculation.
sm=zeros(ng,1);
ksqi=zeros(ng2,1);


      for k = 1:ng2-1
        kdx2 = (pi/ng)*k;
         sm(k) = exp(a1*sin(kdx2)^2 - a2*tan(kdx2)^4);
        ksqi(k) = epsi/((2.0*sin(kdx2)/dx)^2)*sm(k)^2;
     end
       sm(ng2) = exp(a1);
       if a2~=0
          sm(ng2) = 0;
       end
        ksqi(ng2) = +epsi/((2.0/dx)^2)*sm(ng2)^2;

  %} 
%  Transform charge density.   
     rho(1) = rho(1)+rho(ng+1);
     rho(ng+1) = rho(1); 
    
     rhok = fft(rho(1:ng));  
 
    
        tol = 1;
    for k = 0:ng-1
           if k == 0
              ii = tol;
           elseif k <= (ng)/2 && k >= 0
                ii = k;
            elseif k > (ng)/2
                ii = k - (ng);
            end            
            phik(k+1) = -rhok(k+1)/(2*pi*ii/L)^2;

esestot(t+1)=(esestot(t+1)-(phik(k+1)*rhok(k+1)))/2;

    end
    
    phi(1:ng) =-ifft(phik(1:ng));
    phi(ng+1) = phi(1);
    phi=real(phi);

for k=1:mplot+1    
esem(t+1,k)=esem(t+1,k)-phik(k)'*rhok(k)/2;
end

     
  %{  
     hdx = 0.5*dx;
      for j = 1:ng
       rhok(j) = rhok(j)*hdx;
       scrach(j) = 0;
      end
    
     rhok(1)=0;
%  Calculate phik and field energy.
      eses = 0;
      phik(1) = 0;
      for k = 2:ng2
          
       kk = ng - k;
       phik(k) = ksqi(k)*rhok(k);
       phik(kk) = ksqi(k)*rhok(kk);
       eses = eses + rhok(k)*phik(k) + rhok(kk)*phik(kk);
       
      end
      phik(ng2) = ksqi(ng2)*rhok(ng2);
      ese(t+1) = (2.0*eses + rhok(ng2)*phik(ng2))/(2.0*l);
     
%      Save specified mode energies.  
        li = 1.0/l; 
      for km = 1:mplot
       k = km;
       if k~=0 
       kk = ng - k;
       esem(t+1,km) = (rhok(k)*phik(k) + rhok(kk)*phik(kk))*li;
       if( k==kk ) 
       esem(t+1,km) = 0.25*esem(t+1,km);
       end
       end
       end
   
     for k = 2:ng2
       kk = ng - k;
       rhok(k) = sm(k)*rhok(k);
       rhok(kk) = sm(k)*rhok(kk);
     end
       rhok(ng2) = sm(ng2)*rhok(ng2);
       
  %  Inverse transform phi and smoothed rho.   
       for k = 1:ng
         rhok(k) = rhok(k)*li;
         phi(k) = phik(k)*li;
       end
       
       phi =real(-ifft(real(phik))); 
       phi(ng+1) = phi(1);
       rho(ng+1) = rho(1);
        
eses1(t+1)=-phik(2)'*rhok(2)*(2*dx/sum(N));
eses2(t+1)=-phik(3)'*rhok(3)*(2*dx/sum(N));
eses3(t+1)=-phik(4)'*rhok(4)*(2*dx/sum(N));
eses4(t+1)=-phik(5)'*rhok(5)*(2*dx/sum(N));   
%}
  
%E centered differene across two cells
if iw==1||2
for j=2:ng
    E(j,t+1) = (phi(j-1)-phi(j+1))/(2*dx);
end
E(1,t+1) = (phi(ng)-phi(2))/(2*dx);
E(ng+1,t+1) = E(1,t+1);

 %  Centered difference across 1 cell.
elseif iw==3
     dxi = 1.0/dx;
  for j = 1:ng
  E(j,t+1) = (phi(j) - phi(j+1))*dxi;  
  end
  E(ng+1,t+1) = E(1,t+1);
end

Emax(1,t+1)=max(E(:,t+1).^2);
%}
%poisson
%E(:,t+1)= ex;
%  Electric field has not been renormalized yet.
      ael = 1;
a = E(:,t+1);
e2 =E(:,t+1).^2;
EnergiaP(t+1,1)= 0.5*dx*sum(e2);





