function [Det]= HV2_UMN()
% Model of HV2 at UMN
%
%
%--------------------------------------

% let's use the UCB fridge object, since their CP is going to be hot as
% well
fridge = fUCB();
fridge.T_MC=45e-3; % need to check this number

% TES object
W= MaterialProperties('W');

W.Tc=80e-3;%[mK]
Det =HVZIP(fUCB(),eCDMSII(fUCB),'Ge',W,[]);

%Measured Quantities
% NEED TO CHECK ALL OF THESE NUMBERS!
Det.Meas.fPt       = 0.15;
Det.Meas.sigPt0    = 400;%[eV]
Det.Meas.tau_Pabsb = 200e-6;%[us]
