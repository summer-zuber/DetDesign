function fitG115=G115_fit2(lgc_svplt)
% G115_fit2
% Let's fit TES parameters(fOp,gPep_v,Tc,w_Tc) to measured dIdV and IV quantities(loop gain, Tau0, P0,R0,beta)

%Measured Quantities
measPo=2e-12;%[W]
measLG= 6.5;
measbeta=1;
measRo = 85e-3;
meastau0 = 200e-6;%[s]

meas_FileNm='/Users/mpyle1/Dropbox/LAB_Old/R482/Low Bias PSDs/PSD_Transition.mat';
if exist(meas_FileNm,'file')
    lgc_pltmeas=true;
    load(meas_FileNm)
else
    lgc_pltmeas=false;
    display('please supply correct measurement filenm')
end

if nargin < 1 
  lgc_svplt = true;
end    

%1) let's find the fractional operating point
fitG115= HV4mm_4ch;
fitfOp= measRo/fitG115.TES.Rn;

%2) Let's vary 
%       - a scaling factor for gPep_v (the electron phonon coupling) until we find Tau0 [coupled with other 2 equations]
%       - Tc until we find P0 [coupled with other 2 equations]
%       - w_Tc until we find loopgain

guessX=[1,45e-3,5e-4]; %[guess fgPep_v,guess Tc, guess wTc];

options=optimset('fsolve');
options=optimset(options,'MaxFunEvals',2000);

vary= @(x)vary_gPep_Tc_wTc(x,fitfOp,meastau0,measLG,measPo);

fitX=fsolve(vary,guessX,options);
[dY,fitG115]=vary(fitX);

[sigPt_OF,fitG115] = SimulatedNoise_1TES(fitG115);

%let's calculate the DC parasitic power required (the equivalent of a Tbath
%shift)
Tbath_actual = 40e-3;%[k]
fitG115.TES.PpEst= fitG115.TES.gPep_v.* fitG115.TES.vol.*(fitG115.fridge.T_MC.^fitG115.TES.nPep - Tbath_actual.^fitG115.TES.nPep); %[W]


figure(11)
title({'G115 Noise','fitting \Sigma, T_{c},w_{Tc} to dIdV'})
if lgc_pltmeas
    hold on
    for jj=1:4
        plot(Meas.freq,Meas.Si,'-m');
    end
    hold off
    xlim([1e1,1e5])
    ylim([1e0,1e2])
end

if lgc_svplt
    saveas(11,'G115_fit2.png','png')
end
 
end

function [dY,Det] = vary_gPep_Tc_wTc(X,fitfOp,meastau0,measLG,measPo)
    W=MaterialProperties('W');
    W.gPep_v = W.gPep_v.*X(1);
    W.Tc=X(2);
    W.wTc=X(3);
    
    fridge= fUCB;
    
    Det = HV4mm_4ch(fridge,[],[],W,[],fitfOp,[]);
    Det= SimpleEquilibrium_1TES(Det);
    
    dY = [Det.TES.tau0-meastau0; Det.TES.LG-measLG ; Det.TES.Po-measPo];
end

% function dtau0 = Tau0_varyTc(Tc)
%     W=MaterialProperties('W');
%     W.Tc=Tc;
%     
%     Det = HV4mm_4ch([],[],[],W,[],fitfOp,[]);
%     Det= SimpleEquilibrium_1TES(Det);
%     
%     dtau0 = Det.TES.tau0-meastau0;
% end