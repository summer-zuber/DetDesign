ltrap = 180e-6;
ltes  =  45e-6;
fwal  =  0.5;

lwal= 2*ltes*fwal;
lfin = .75*ltrap;

eff_absb=Effqp_Ideal_cfin(lfin,lwal,ltrap,true);


nfin=100;
lfin=linspace(.01,3,nfin)'*ltrap;
x=lfin./ltrap;

Amax = 4*lfin(end).^2;

%let's calculate Rin and Rout
Rin  = lwal/(2*pi);
Rout = Rin + lfin;
Afin= pi*(Rout.^2-Rin.^2);

%

eff_absb=Effqp_Ideal_cfin(lfin,lwal,ltrap,false);


figure(11)
plot(x,eff_absb*100,'-k')
hold on
plot(x, eff_absb(end).*(Afin(end)./Afin)*100,'--b')
hold off
title('QP collection efficency for 2D fins')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} %')
ylim([0,100])
legend({'Estimate','1/A_{Al}'})

figure(12)
plot(x,(eff_absb.*Afin)/Amax,'-k')
title('First Phonon Wave Collection Capability for 2D fins')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} A_{fin} ')
%ylim([0,1])