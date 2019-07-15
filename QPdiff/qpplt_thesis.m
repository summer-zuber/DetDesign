ltrap=180e-6;
lfin=[0:.01:4]'*ltrap;

[eff_absb]=Effqp_absb_rfin(lfin,ltrap);
%[eff_absb2]=Effqp_absb_rfin_old(lfin,ltrap*sqrt(2));

x= lfin./ltrap;

figure(11)
plot(x,eff_absb*100,'-k')
hold on
%plot(lfin./ltrap,eff_absb2*100,':g')
plot(x, 1./x*100,'--b')
hold off
title('QP collection efficency for 1D simple fins')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} %')
ylim([0,100])

figure(12)
plot(x,eff_absb.*x,'-k')
title('First Phonon Wave Collection Capability')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} l_{fin}/l_{trap} ')
ylim([0,1])