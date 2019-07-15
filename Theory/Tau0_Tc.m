%Tau0_Tc

%let's make a plot of Tau0 as a function of Tc for our system:
nTc=1e3;
Tc= linspace(20e-3,100e-3,nTc)';

W=MaterialProperties('W');

C=W.fCsn .* W.gC_v .* Tc;%[J/K/m^3]
G=W.nPep .* W.gPep_v .*Tc.^(W.nPep-1);%[W/K/m^3]
tau0=C./G;%[s]

%G115 measured tau0:
tau0_G115 = [215e-6,212e-6,197e-6,197e-6,219e-6,203e-6,213e-6,201e-6];%[s]

figure(1)
plot(Tc*1e3,tau0*1e6,'-k')
hold on
for j=1:8
    plot(Tc*1e3,tau0_G115(j)*1e6*ones(nTc,1),'--b')
end
hold off
xlabel('T_{c}[mK]')
ylabel('Thermal Fall Time (C/G) [us]')
xlim([20,100])
ylim([100,2000])
set(gca,'xscale','log','yscale','log')
grid on
hleg=legend({'Expected','G115 Fit'},'location','best')
set(hleg,'FontSize',24)
title('Thermal Fall Time vs T_{c}')
saveas(1,'tau0_Tc.png','png')