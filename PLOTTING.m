for species=1:nsp
figure
plot(gridt,ke(:,species),gridt,real(EnergiaP),gridt,real(te),'LineWidth',2)
set(gca,'FontSize',12)
legend('kinetic', 'electric potential', 'total')
grid on 
xlim([dt nt*dt])
title([example ,'Energies'],'FontSize',16)
xlabel('Time (arb. units)','FontSize',16)
ylabel('Energy (arb. units)','FontSize',16)
end
saveas(gcf,[sprintf(example), 'Energies.png'])

figure
for k=2:mplot+1
   txt=['Mode ',num2str(k-1),' Energy'];
semilogy(gridt,real(esem(:,k)),'DisplayName',txt,'LineWidth',2)
set(gca,'FontSize',12)
 grid on 
%ylim([1*10^-5 5]);
xlim([3*dt nt*dt]);
%ylim([10^-15 10^-12])
title([example ,'Electric Mode Energies with Excitation k= ' num2str(mode(1))])
%subplot(3,1,3);
hold on
xlabel('Time (arb. units)','FontSize',16)
ylabel('Energy (arb. units)','FontSize',16)
end
hold off
legend show
saveas(gcf,[sprintf(example), 'Mode Energies with Excitation k= ' num2str(mode(1)),'.png'])

figure
semilogy(gridt,EnergiaP,'LineWidth',2)
 grid on 
xlim([dt nt*dt]);
xlabel('Time (arb. units)','FontSize',16)
ylabel('Energy (arb. units)','FontSize',16)
title([example ,'Field Energy'],'FontSize',16)
saveas(gcf,[sprintf(example), 'Field Energy.png'])
figure

for species=1:nsp
subplot(nsp,1,species)
plot(gridt,ke(:,species),'LineWidth',2)
set(gca,'FontSize',12)
grid on
title([example ,'KE Energy for Species v0  = ', num2str(v0(species))],'FontSize',16)

xlim([0 nt*dt]);
xlabel('Time (arb. units)','FontSize',16)
ylabel('Energy (arb. units)','FontSize',16)
hold off
end
saveas(gcf,[sprintf(example), 'KE.png'])


figure
for species=1:nsp
subplot(nsp,1,species)
plot(gridt,de(:,species),'LineWidth',2)
set(gca,'FontSize',12)
grid on
title([example ,'Drift Energy for Species v0 = ', num2str(v0(species))],'FontSize',16)
xlim([0 nt*dt]);

xlabel('Time (arb. units)','FontSize',16)
ylabel('Energy (arb. units)','FontSize',16)
hold off
end
saveas(gcf,[sprintf(example), 'DE.png'])


figure
for species=1:nsp
subplot(nsp,1,species)
plot(gridt,therme(:,species),'LineWidth',2)
set(gca,'FontSize',12)
grid on
title([example ,'Thermal Energy for Species v0 = ', num2str(v0 (species))],'FontSize',16)
xlim([0 nt*dt]);

xlabel('Time (arb. units)','FontSize',16)
ylabel('Energy (arb. units)','FontSize',16)
hold off
end
saveas(gcf,[sprintf(example), 'TE.png'])




