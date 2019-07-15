function [Det] = Heating_Estimates(Det);
% Let's calculate the heat load due to the Detector.
%
% IMPORTANT: WE EXPLICITLY ASSUME THAT ALL RP is on the TES leg of the
% circuit
%
%
% 12/9/13 MCP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\\\\\\\\\\ Calculate Bias Temperature \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Det = TES_Dynamics_Simple(Det);

%\\\\\\\\\\ Heating Estimates \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%1) TES Heating 
Heat=[];

Heat.TES = Det.TES.gPep_v .* Det.TES.vol .* (Det.TES.Tbias.^Det.TES.nPep - Det.fridge.T_MC.^Det.TES.nPep); %[W]

%2) Rp Heating
%We need to calculate the TES current:
Ites= sqrt(Heat.TES/Det.TES.Ro); %[A]

Heat.Rp = Ites.^2 * Det.elec.Rp; %[W]

%3) Rshunt Heating
%We need to calculate the shunt current:
Ishunt = Ites .* (Det.elec.Rpt+Det.TES.Ro)/Det.elec.Rs; %[A]

Heat.Rs = Ishunt.^2 .*Det.elec.Rs; %[W]

%\\\\\\\\\\\ Find Total Heat Loads for Detector \\\\\\\\\\\\\\\\\\\\\\\\\\\
Heat.TES = Det.nP * Heat.TES;%[W]
Heat.Rp  = Det.nP * Heat.Rp; %[W]
Heat.Rs  = Det.nP * Heat.Rs; %[W]

%\\\\\\\\\\\ Add to Det \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Det.Heat=Heat;

