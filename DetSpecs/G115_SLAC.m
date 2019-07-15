function G115=G115_SLAC(lgc_svplt)
%G115_SLAC

%Measured Quantities
measPo   = 5e-12;%[W]
measLG   = 6.5*2.5;
measbeta = 0.0;
measRo   = 190e-3*0.50;%[Ohm]
meastau0 = 200e-6;%[s]
measRn   = 190e-3; %[Ohm]

%1) let's find the fractional operating point
W        = MaterialProperties('W');
W.Tc     = 45e-3;%[K]
W.wTc    = 2.5e-4;%[K]

%let's bump up the electron-phonon coupling by x2.5
W.gPep_v = W.gPep_v*2.5;


G115     = HV4mm_4ch(fSNOLAB,eSLAC_G115,'Ge',W,[],measRo/measRn);
G115.TES.beta= measbeta;

[sigPt_OF,G115] = SimulatedNoise_1TES(G115);

figure(11)
xlim([1e1,1e6])
title('50% Rn G115 @SLAC \beta=0')