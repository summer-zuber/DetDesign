function [ResPt,Det] = ResPt_Estimate(Det);
% Let's estimate the OF total energy resolution for devices with these
% assumptions:
%
%   1) alpha/n >> 1
%   2) L/Ro >> tau_etf & tau_phonon
%   3) Squid Noise is negligible
%
% Outputs:
%
%   1) theoretical prediction for sigPt
%

%%%%%%%%% Physical Constants and Invariant Characteristics %%%%%%%%%%%%%%%%
pc=PhysicalConstants();

%%%%%%%%%%%%%%%%%%%% TES Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Det=SimpleEquilibrium_1TES(Det);

%%%%%%%%%%%%%%%%%%% TFN Noise Coefficient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%let's use the ballistic version
fNtfn= ((Det.fridge.T_MC/Det.TES.To).^(Det.TES.nPep+1) +1)/2;

%let's save the TFN noise for a table
Det.TES.Sp_tfn = 4*pc.kb .* Det.TES.To.^2.*Det.TES.Gep .* fNtfn;

%%%%%%%%%%%%%%%%%%% Energy Resolution Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%
%------ Theoretical Resolution Calculation --------------------------------
% Shunt Resistor and Parasitic Resistor Noise is very important let's
% put into the calculation

% let's calculate an effective temperature
Tstar = (Det.elec.Tl.*Det.elec.Rl).*Det.TES.Ro./(Det.TES.Ro-Det.elec.Rl).^2;

% let's calculate an effective sensor bandwidth
w_star = sqrt((fNtfn.*Det.TES.nPep.*Det.TES.Tc + Tstar)./(Det.TES.Tc+Tstar)).*Det.TES.w_etf;

ResPt = sqrt(Det.nP * 4*pc.kb.*Det.TES.Tc.^2.*Det.TES.Gep./Det.eEabsb.^2 .* ((fNtfn.*Det.TES.nPep.*Det.TES.Tc + Tstar)./(Det.TES.nPep.*Det.TES.Tc)) .*(w_star+Det.w_Pabsb)./w_star./Det.w_Pabsb)/pc.qe;
Det.ResPt=ResPt;

% Let's assume that the Rs and Rp johnson noise doesn't
% matter
Det.ResPt_simple= sqrt(Det.nP.* 4.*pc.kb.*Det.TES.To.^2.*Det.TES.Gep./Det.eEabsb.^2 .*fNtfn.*(sqrt(fNtfn.*Det.TES.nPep)*Det.TES.w_etf+Det.w_Pabsb)./(sqrt(fNtfn.*Det.TES.nPep).*Det.TES.w_etf)./Det.w_Pabsb)/pc.qe;

