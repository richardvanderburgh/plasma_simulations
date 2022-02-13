%ES1 - Main Program
clear variables
close all
clc

FIRST_EE
INIT
t=0;
SETRHO
FIELDS
SETV
    
axis tight manual % this ensures that getframe() returns a consistent size
v = VideoWriter(sprintf(example));
open(v);

for t=1:nt
   
    if mod(t,ixvx)==0|| t==1 && ixvx~=0  
     if t==1 
     phasecounter=1; 
     PhaseSpacePlots=figure;
     hold on
     end
    for species=1:nsp
    set(0,'CurrentFigure',PhaseSpacePlots)
    
    %h1=subplot(3,3,phasecounter) ;
    
    %subplot(2,2,1)
    %scatter(x(:,species)/dxi,vx(:,species)*dxdt, '.'); 
    scatter(x(:,species)/dxi,vx(:,species), '.'); 
    %set(gca,'FontSize',12)
    grid on 
    %ylim([-x1(1), x1(1)]);    
    ylim([-3, 3]);    
    xlim([0 L]);
    hold on
    end
    hold on
    if t==1
    title([sprintf(example), 'Phase Space at Time = ',num2str((t-1)*dt-dt/2),]);
    xlabel('Position (arb. units)');
    ylabel('Velocity (arb. units)');
    elseif t ~= 1 
    title([sprintf(example), 'Phase Space at Time = ',num2str((t)*dt-dt/2),]);
    xlabel('Position (arb. units)');
    ylabel('Velocity (arb. units)');
    end
    phasecounter=phasecounter+1;
    hold off
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
    end
    if t==nt && ixvx~=0 
    set(0,'CurrentFigure',PhaseSpacePlots)
    %suptitle([sprintf(example), 'Phase Space'])
    saveas(PhaseSpacePlots,[sprintf(example), 'PhaseSpace.png'])
    end
   
    ACCEL
    MOVE
    FIELDS
    
    for species=1:nsp
    te(t,1)=te(t,1)+ke(t,species);
    end
    te(t,1)=te(t,1)+EnergiaP(t,1);
 
   if mod(t,irho)==0|| t==1 && irho~=0;
   
       
        if t==1
    rhocounter=1;  
    ChargeDensityPlots=figure;
    hold on
        end
    set(0,'CurrentFigure',ChargeDensityPlots)
    h2=subplot(3,3,rhocounter) ;
     %subplot(2,2,2)
    plot(gridx,real(rho),'LineWidth',2)
    xlim([0 L]);
    ylim([-4 -0]);
    grid on 
    hold on
    xlabel('Position');
    ylabel('Charge Density');
    if t ==1
    title(['Time = ',num2str((t-1)*dt)]);  
    elseif t ~= 1 
    title(['Time = ',num2str((t)*dt)]);
    xlabel('Position (arb. units)');
    ylabel('Charge Density (arb. units)');
    grid on 
    end
    rhocounter=rhocounter+1;
    drawnow
    %hold off
    end
    if t==nt && irho~=0
    set(0,'CurrentFigure',ChargeDensityPlots)
    suptitle([sprintf(example), 'Charge Density'])
    saveas(ChargeDensityPlots,[sprintf(example), 'ChargeDensity.png'])
    end
   
   if mod(t,iphi)==0 || t==1 && iphi~=0 
   
       if t==1
       phicounter=1; 
       PotentialPlots=figure;
       hold on
       end
        set(0,'CurrentFigure',PotentialPlots)
    h3=subplot(3,3,phicounter) ;
   
    %subplot(2,2,3)
    plot(gridx,real(phi),'LineWidth',2)
    xlim([0 L]);
    %ylim([-3 -1]);
     grid on 
    %hold on
    xlabel('Position (arb. units)');
    ylabel('Potential (1arb. units)');
    if t ==1
    title(['Time = ',num2str((t-1)*dt)]);
    elseif t ~= 1
    title(['Time = ',num2str((t)*dt)]); 
    xlabel('Position (arb. units)');
    ylabel('Potential (arb. units)');
    grid on
    end 
    drawnow
    phicounter=phicounter+1; 

    hold off
   end
    if t==nt && iphi~=0
    set(0,'CurrentFigure',PotentialPlots)
    suptitle([sprintf(example), 'Electric Potential'])
    saveas(PotentialPlots,[sprintf(example), 'ElectricPotential.png'])
    end
   
  if mod(t,iE)==0 || t==1 && iE~=0 
    
   if t==1
       Ecounter=1; 
       ElectricFieldPlots=figure;
       hold on
   end
    set(0,'CurrentFigure',ElectricFieldPlots)
    h4=subplot(3,3,Ecounter) ;
        
    %subplot(2,2,4)
    plot(gridx,real(E(t,:)),'LineWidth',2)
    xlim([0 L]);
    ylim([-1 1]);
     grid on 
    %hold on
    xlabel('Position (arb. units)');
    ylabel('Electric Field (1arb. units)');
    if t ==1
    title(['Time = ',num2str((t-1)*dt)]);
    elseif t ~= 1
    title(['Time = ',num2str((t)*dt)]); 
    xlabel('Position (arb. units)');
    ylabel('Electric Field (arb. units)');
    grid on
    end 
    Ecounter=Ecounter+1;
    drawnow
    hold off
    
  end
   
   if t==nt && iE~=0
  set(0,'CurrentFigure',ElectricFieldPlots)
  suptitle([sprintf(example), 'Electric Field'])
  saveas(ElectricFieldPlots,[sprintf(example), 'ElectricField.png'])
   end

    if mod(t,ifvx)==0 || t==1 && ifvx~=0 
      
    if t==1
    fvxcounter=1;  
    fvxPlots=figure;
    hold on
    end
    
    set(0,'CurrentFigure',fvxPlots)
    h5=subplot(3,3,fvxcounter) ;
  
    hist(vx);
    grid on
   
    if t==1
    title(['fvx at time ',num2str((t-1)*dt-dt/2),]);
    xlabel('velocity');
    
    elseif t ~= 1 
    title(['fvx at time ',num2str((t)*dt-dt/2),]);
    xlabel('velocity');
    %xlim([-3 3]);
    %ylim([0 50]);
    end
    fvxcounter=fvxcounter+1;
    drawnow
    hold off
    end 
   
  if t==nt && ifvx~=0
  set(0,'CurrentFigure',fvxPlots)
  suptitle([sprintf(example), 'Velocity Distribution'])
  saveas(fvxPlots,[sprintf(example), 'Velocity Distribution.png'])
   end
  
    %t*dt-dt/2
   
end
close(v);
%plotspectrnew
wk_bland
%wk_zpic


PLOTTING

%End of Simulation :)