function [Det]= S12C_R466()
% Model of S12C in Run 466
%
% Run 466: 
%          -in R466 S12C had all channels operating
%          - side 1 already low Tc
%          - side 2 ion-implanted
%          - Tc ~ 65mK ->  NEED TO CHECK .XLSX sheets ... is there a
%          sizable difference in Tc between the 2 sides? Also, note that in
%          R472 we measured much higher Tc's confusing.
%
%--------------------------------------

% fridge object
f75uW_R466 = fUCB();
f75uW_R466.T_MC=36e-3; % need to check this number

% TES object
W= MaterialProperties('W');

%W.Tc=65e-3;%[mK] -> this was measured in R466
W.Tc=72e-3;%[mK] -> this was measured in R472
Det =iZIP4(fUCB(),eCDMSII(fUCB),'Si',W,[]);

%We only used 7 channels rather than 8
Det.nP=7;

%Measured Quantities
Det.Meas.fPt       = 0.1925;
Det.Meas.sigPt0    = 48.6;%[eV]
Det.Meas.tau_Pabsb = 175e-6;%[us]
