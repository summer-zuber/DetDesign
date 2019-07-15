function fitG115=G115_2DOFfit1(lgc_svplt)
%G115_fit1

%Initialize
fitG115 = HV4mm_4ch([],[],[],[],[],[],[],true,[]);
fitG115 = SimpleEquilibrium_2TESParallel(fitG115);

% //////// Measured Quantities ///////////
%Let's try to fit channel B [red] exactly (http://titus.stanford.edu/cdms_restricted/detector_physics/HV/ebook/170511c/)

meas=[];

%Measured Bath Quantities
Meas.T_MC= 40e-3;%[K]

%Measured Po(T_{MC}) Quantities
meas.gPep_v= 0.37e9; %[W/m^3/K^5] ->  from G23R
meas.nPep=5;

%Measured IV Quantities
meas.Ro = 112e-3;%[Ohm] -> used
meas.Po = 2e-12;%[W]

%Measured dIdV Quantities
meas.A = 145e-3; %[Ohm] %- > used to find beta
meas.B= -0.4; %-> used to find LG
meas.C= - 0.15; 
meas.tau_el = 0.6e-6;%[us] %-> used to find Ltot

meas.tau_I  = -100e-6;%[s]

meas.tau_1  =  500e-6;%[s]

meas.C = .15; %- used to calculate Fg

%Measured dIdV Quantities
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


%------ Basic Fit Parameters -------------------

fit=[];

%1) let's find the fractional operating point
fit.fOp  = meas.Ro/fitG115.TES.Rn;

%2) Lets find other parameters:
fit.beta = (meas.A-meas.Rl)./meas.Ro-1;

%3) Let's find the effective loop gain:
fit.LG = meas.B ./(meas.B+meas.Ro.*(2+fit.beta));

%4) Let's measure the total inductance
fit.Lt = meas.tau_el .* meas.A; %[H]

%5) Thermal Conductance ratios which change the magnitude of the 2nd pole
%expression
fit.Fg = fit.C .* (1-Meas.LG);

%-----  Fitting ---------------------------------
% Equations:
% 1) tau_I  
% 2) tau_1
% 3) Po(Tb)
% 4) dPo/dTb
% 5) F


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

function [dY,fitG115] = vary_Tc_wTc_Qp(X,fitfOp,meas)
    % X(1) = TES Tc [K]
    % X(2) = TES transition width [K]
    % X(3) = Parasitic Power [pW]

    W=MaterialProperties('W');
    W.Tc=X(1);
    W.wTc=X(2);
    
    fridge = fUCB(meas.Tbath);
    
    fitG115 = HV4mm_4ch(fridge,[],[],W,[],fitfOp,X(3)*1e-12,[]);
    fitG115= SimpleEquilibrium_1TES(fitG115,meas.beta);
    
    dY = [(fitG115.TES.tau0-meas.tau0)*1e6; fitG115.TES.LG-meas.LG ; (fitG115.TES.Po-meas.Po)*1e12];
end