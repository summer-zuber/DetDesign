ltes  = 150e-6;
ltrap = 180e-6;
lfin  = [.1:.01:4]'*ltrap;

fwal  = 1/400; 
lscat = 150e-9;

%Derived Quantities
Rin= ltes./(2*pi);
Rout= Rin+lfin;
Afin_C= pi*(Rout.^2-Rin.^2);

Afin_R = 2*ltes.*lfin;

[effabsb_R]=Effqp_absb_rfin(lfin,ltrap);
[effabsb_R2]=Effqp_absb_rfin(lfin,ltrap,fwal,lscat);
[effabsb_C]=Effqp_absb_cfin(lfin,ltes,ltrap,fwal,lscat);

x= lfin./ltrap;

figure(11)
plot(x,effabsb_R*100,'-k')
hold on
plot(x,effabsb_R2*100,'-g')
plot(x,effabsb_C*100,'-m')
plot(x,25*mean(Afin_R)./Afin_R,':c')
hold off
title('QP collection efficiency')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} %')
legend({'1D ideal','1D f=1/400','2D f=1/400'})
ylim([0,100])

%let's multiply by the total area of the fin:

figure(12)
plot(x,effabsb_R.*Afin_R*1e12,'-k')
hold on
plot(x,effabsb_R2.*Afin_R*1e12,'-g')
plot(x,effabsb_C.*Afin_C*1e12,'-m')
hold off
title('First Phonon Wave Collection Capability')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} A_{fin} (um^2)')
legend({'1D ideal','1D f=1/400','2D f=1/400'})
