function fitG115 = G115_fit4(lgc_svplt)
%G115_fit4

%Measured Quantities
meas =[];
meas.Po=2e-12;%[W]
meas.LG= 6.5;
meas.beta=1;
meas.Ro = 85e-3;
meas.tau0 = 200e-6;%[s]
meas.Tbath = 40e-3;%[K]

meas_FileNm='/Users/mpyle1/Dropbox/LAB_Sadoulet/R482/Low Bias PSDs/PSD_Transition.mat';
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

fit=[];
%1) let's find the fractional operating point
fitG115 = HV4mm_4ch;
fitfOp  = meas.Ro/fitG115.TES.Rn;

%2) Let's vary 
%       - Tc until we find Tau0
%       - w_Tc until we find loopgain
%       - Qp until we find the correct P0

%Notice that the parasitic power is in units of pW ... if not the fitter
%doesn't work :(
guessX=[80e-3,5e-4,10]; %[guessTc,guesswTc,guessQp];


options=optimset('fsolve');
options=optimset(options,'MaxFunEvals',5000,'MaxIter',5000);

vary= @(x)vary_Tc_wTc_Qp(x,fitfOp,meas);

fitX=fsolve(vary,guessX,options);
[dY,fitG115]=vary(fitX);

[sigPt_OF,fitG115] = SimulatedNoise_1TES(fitG115);

figure(11)
title({'G115 Noise','fitting T_{c},T_{bath},w_{Tc} to dIdV'})
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
    saveas(11,'G115_fit1.png','png')
end
 
end

function [dY,Det] = vary_Tc_wTc_Qp(X,fitfOp,meas)
    % X(1) = TES Tc [K]
    % X(2) = TES transition width [K]
    % X(3) = Parasitic Power [pW]

    W=MaterialProperties('W');
    W.Tc=X(1);
    W.wTc=X(2);
    
    fridge = fUCB(meas.Tbath);
    
    Det = HV4mm_4ch(fridge,[],[],W,[],fitfOp,X(3)*1e-12,[]);
    Det= SimpleEquilibrium_1TES(Det,meas.beta);
    
    dY = [(Det.TES.tau0-meas.tau0)*1e6; Det.TES.LG-meas.LG ; (Det.TES.Po-meas.Po)*1e12];
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