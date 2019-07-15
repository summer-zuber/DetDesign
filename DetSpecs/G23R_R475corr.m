function [Det]= G23R_R475corr()
% Model of G23R in R475 after correction for aSi removal
%
% G23R is a 1 sided iZIP. The other side has no Al. However it has 100% aSi
% coverage. From looking at the athermal phonon fall times, we know that
% there is a big athermal phonon sink. We think it could be from the aSi.
% Here we correct the numbers for this.
%--------------------------------------

% fridge object
f75uW_R466 = fUCB();
f75uW_R466.T_MC = 36e-3; % need to check this number

% TES object
W= MaterialProperties('W');

W.Tc = 52e-3;%[mK]
Det  = iZIP4_1sided(fUCB(),eCDMSII(fUCB),'Ge',W,[]);

%let's correct all the measured quantities to account for the fact that the
%aSi sunk a bunch of athermal phonons
tau_PabsbExp= 775e-6*2;%[s]
w_PabsbExp= 1./tau_PabsbExp; %[1/s]

effaSi= w_PabsbExp/Det.w_Pabsb;

%"Corrected" Quantities
Det.Meas.fPt       = 0.112 / effaSi;
Det.Meas.sigPt0    = 52.*sqrt(effaSi);%[eV]
Det.Meas.tau_Pabsb = tau_PabsbExp;%[s]
Det.Meas.w_Pabsb   = w_PabsbExp;%[1/s]  
