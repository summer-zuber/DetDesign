function [Det]= G23R_R475()
% Model of G23R in R475
%
% G23R is a 1 sided iZIP. The other side has no Al ... it's just aSi
%
%
%
%--------------------------------------

% fridge object
f75uW_R466 = fUCB();
f75uW_R466.T_MC = 36e-3; % need to check this number

% TES object
W= MaterialProperties('W');

W.Tc = 52e-3;%[mK]
Det  = iZIP4_1sided(fUCB(),eCDMSII(fUCB),'Ge',W,[]);

%Measured Quantities
Det.Meas.fPt       = 0.112;
Det.Meas.sigPt0    = 52;%[eV]
Det.Meas.tau_Pabsb = 875e-6;%[s]
Det.Meas.w_Pabsb   = 1/Det.Meas.tau_Pabsb;%[1/s]
