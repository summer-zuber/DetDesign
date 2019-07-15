%G115f3 = G115_fit3();
fitdIdV = mk_fake_fitdIdV(G115f3);

%let's vary Ro by +/-100mOhms
nvary=100;

Ro= fitdIdV.Ro + linspace(-10e-3,10e-3,nvary)';

LG=zeros(nvary,1);
tau0=zeros(nvary,1);
for jj=1:nvary
    fitdIdV.Ro = Ro(jj);
    TES = estTES_fitdIdV(fitdIdV);
    LG(jj)=TES.LG;
    tau0(jj)=TES.tau0;
end

figure(1)
plot(Ro*1e3,tau0*1e6,'-k')
xlabel('R_{0} [mOhms]')
ylabel('tau_{0} [us]')
title('Variation of C/G with R_{0} systematics')
saveas(1,'varyTau0_Ro.png','png')

figure(2)
plot(Ro*1e3,LG,'-k')
xlabel('R_{0} [mOhms]')
ylabel('tau_{0} [us]')
title('Variation of C/G with R_{0} systematics')
saveas(2,'varyLG_Ro.png','png')

%let's put Ro back to normal
fitdIdV.Ro=G115f3.TES.Ro;

%--------------------------------------------------------------------------

%Let's vary dIdV by 10mOhm

dVdI0= 1./fitdIdV.dIdV0 + linspace(-10e-3,10e-3,nvary)';

LG=zeros(nvary,1);
tau0=zeros(nvary,1);
for jj=1:nvary
    fitdIdV.dIdV0 = 1./dVdI0(jj);
    TES = estTES_fitdIdV(fitdIdV);
    LG(jj)=TES.LG;
    tau0(jj)=TES.tau0;
end

figure(11)
plot(dVdI0*1e3,tau0*1e6,'-k')
xlabel('dVdI(0) [mOhms]')
ylabel('tau_{0} [us]')
title('Variation of C/G with dVdI(0) systematics')
saveas(11,'varytau0_dIdV0.png','png')

figure(12)
plot(dVdI0*1e3,LG,'-k')
xlabel('dVdI(0) [mOhms]')
ylabel('tau_{0} [us]')
title('Variation of C/G with dVdI(0) systematics')
saveas(12,'varytau0_dIdV0.png','png')

%---------------------------------
