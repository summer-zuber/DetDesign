ltrap = 180e-6*sqrt(2);
ltes  = 100e-6;
fwal  = 0.5;

nfin=200;

lfin2=linspace(.01,3,nfin)'*ltrap;
lfin1=linspace(.01,7,nfin)'*ltrap;

x_2D=lfin2./ltrap;
x_1D=lfin1./ltrap;

Amax = pi*lfin2(end).^2;

Alim = pi*(2*ltrap)^2;

%-------------- Rectangular Fin -------------------------------------------
Afin_1D= 2.*ltes.*lfin1;

effabsb_1D=Effqp_Ideal_rfin(lfin1,ltrap);

%-------------- Circular Fin ----------------------------------------------
lwal= 2*ltes*fwal;

%let's calculate Rin and Rout
Rin  = lwal/(2*pi);
Rout = Rin + lfin2;
Afin_2D= 2*ltes.*lfin2+ pi*lfin2.^2;
Afin_2Dchk= pi*(Rout.^2-Rin.^2);


effabsb_2D=Effqp_Ideal_cfin(lfin2,lwal,ltrap,false);

%------------ Plots -------------------------------------------------------
figure(10)
plot(Afin_1D*1e12,effabsb_1D*100,'-k')
hold on
plot(Afin_2D*1e12,effabsb_2D*100,'-b')

plot(Afin_1D*1e12, effabsb_1D(end).*(Afin_1D(end)./Afin_1D)*100,':k')
plot(Afin_2D*1e12, effabsb_2D(end).*(Afin_2D(end)./Afin_2D)*100,':b')
hold off
title('QP collection efficency for 1D & 2D fins')
xlabel('A_{fin} (\mum^2)')
ylabel('\epsilon_{collect} %')
xlim([0,e5])
ylim([0,100])
legend({'1D','2D','1D: 1/A_{fin}','2D: 1/A_{fin}'},'location','northeast')
saveas(10,'plt/QPeff_Afin_1D2D.pdf','pdf')

figure(11)
plot(x_1D,effabsb_1D*100,'-k')
hold on
plot(x_1D, effabsb_1D(end).*(Afin_1D(end)./Afin_1D)*100,':k')

plot(x_2D,effabsb_2D*100,'-b')
plot(x_2D, effabsb_2D(end).*(Afin_2D(end)./Afin_2D)*100,':b')
hold off
xlim([0,3.5])
title('QP collection efficency for 1D & 2D fins')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} %')
ylim([0,100])
legend({'1D','1D: 1/A_{fin}','2D','2D: 1/A_{fin}'},'location','northeast')
saveas(11,'plt/QPeff_lfin_1D2D.pdf','pdf')


figure(12)
plot(x_1D,(effabsb_1D.*Afin_1D)/Amax,'-k')
hold on
plot(x_2D,(effabsb_2D.*Afin_2D)/Amax,'-b')
hold off
xlim([0,3.5])
title('First Phonon Wave Collection Capability for 1D & 2D fins')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} A_{fin} ')
legend({'1D','2D'},'location','east')
saveas(12,'plt/QPeffAfin_lfin_1D2D.pdf','pdf')


figure(13)
plot(x_1D, effabsb_1D.*sqrt(Afin_1D/Amax),'-k')
hold on
plot(x_2D, effabsb_2D.*sqrt(Afin_2D/Amax),'-b')
%plot(x_2D, effabsb_2D.*sqrt(Afin_2Dchk/Amax),':b')
hold off
title('Optimum Filter Resolution for 1D & 2D fins')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect}A_{fin}^{1/2} ')
legend({'1D','2D'},'location','northeast')
saveas(13,'plt/OFres_lfin_1D2D.pdf','pdf')


figure(14)
plot(effabsb_1D.*sqrt(Afin_1D/Amax), effabsb_1D.*Afin_1D/Amax,'-k')
hold on
plot(effabsb_2D.*sqrt(Afin_2D/Amax), effabsb_2D.*Afin_2D/Amax,'-b')
hold off
title('Position and Energy Metrics')
xlabel('Energy: \epsilon_{collect} A_{fin}^{1/2}')
ylabel('Postion: \epsilon_{collect} A_{fin}')
legend({'1D','2D'},'location','northwest')
saveas(14,'plt/OFres_Posres_1D2D.pdf','pdf')
