function [Det] = SimpleEquilibrium_2TESParallel(Det,beta,Qp)
% The dIdV seem to show additional fall time poles which seem roughly
% consistent with some additional internal DOF.
%
% The model that we are using to explain these DOF is the "parallel" model,
% in which the TES has 2 heat capacities (the TES and the Fin Connectors), both of which have thermal
% conductances to the bath.
%
% It (and the convention that we are using) can be found in: 
%
% I. J. Massilita,
% "Complex impedance, responsivity and noise of transition-edge sensors: analytical
% solutions for two- and three-block thermal models"
% arXiv: 1205.5693
%
%
%
% A Quick/Inaccurate Estimate of the TES equilibrium point. Basically, we
% calculate the bias point assuming that beta is equal to zero, and then we
% adhoc add a beta after the fact.
%
% It's a bit sketchy since it doesn't fit the functional form that we assume for TES resistance,
% However, it will produce an dynamical system that will exactly match all small signal measurements of a simple 1 block TES
%   
%   Assumptions:
%       1) beta = 0;f
%       2) R(T)= Rn(1+tan((T-Tc)/w_tc))/2
%
%
% Inputs:
%   1) Detector Object
%   2) beta [Default = 0]
%   3) Qp (parasitic power into the TES)[Default = Whatever is in Detector Object]
%
% 12/9/13 MCP
% 17/05/20 MCP: Added DC parasitic heating source.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(beta)
    %beta=1;
    beta=0;
end

if nargin <3 || isempty(Qp)
    %If the detector object is carrying a parasitic heating variable, then
    %we just leave it ... if not let's add a field
    if ~isfield(Det.TES,'Qp')
        Qp=0;
    end
else
  Det.TES.Qp = Qp;%[W]  
end

%\\\\\\\\\\\ Physical Constants and Invariant Characteristics \\\\\\\\\\\\\
pc=PhysicalConstants();

%------- Bias Point Temperature ------------------------------------
% let's calculate the temperature of the operating point resistance.
% [Notice that if the resistance of the TES changes with current, then this doesn't work]

zeta_o= log(Det.TES.fOp/(1-Det.TES.fOp))/2;
Det.TES.Tto = zeta_o * Det.TES.wTc + Det.TES.Tc; %K


%------- Alpha/Beta at Transition Point -----------------------------------
% let's calculate bias location along the normalized transition curve

%let's also calculate alpha 
Det.TES.alpha= 2.*Det.TES.Tto/Det.TES.wTc ./ exp(zeta_o) ./ (exp(zeta_o)+exp(-zeta_o));
Det.TES.beta = beta;

%------ TES Properties @ Equilibrium --------------------------------------
% TES heat capacity
Det.TES.Ct   =  (Det.TES.fCsn*Det.TES.gC_v.* Det.TES.Tto) * Det.TES.volTES; %J/K

% TES phonon/electron thermal coupling to the bath 
Det.TES.Ktb = Det.TES.gPep_v*Det.TES.volTES;
Det.TES.Gtb = Det.TES.nPep.* Det.TES.Ktb .* Det.TES.Tto.^(Det.TES.nPep-1); %W/K

% Internal Thermal Conductance: Evaluated at TES temperature
Det.TES.Gt1_Tt = Det.TES.nt1.*Det.TES.Kt1 .* Det.TES.Tto.^(Det.TES.nt1-1);

%----- Fin Connector Properties (Additional Heat Capacity: C1) -----------------------
% First we need to calculate, T1, the effective temperature of the fin
% connector. We do this by having the Heat flow into the fin connector
% equal to the heat flow out of the fin connector:

Det.TES.K1b = Det.TES.gPep_v*Det.TES.veff_FinCon*Det.TES.volFinCon;

dE1dt = @(T) (Det.TES.Kt1.*(Det.TES.Tto.^Det.TES.nt1-T.^Det.TES.nt1) - Det.TES.K1b.*(T^Det.TES.nPep-Det.fridge.T_MC.^Det.TES.nPep))*1e12; %[pW]

%let's make a voltage divider out of the thermal impedances as a first
%guess for T1
G1b_Tb = Det.TES.nPep.*Det.TES.K1b.*Det.fridge.T_MC.^(Det.TES.nPep-1);

T1guess = Det.fridge.T_MC + (Det.TES.Tto - Det.fridge.T_MC).* (1./G1b_Tb)./(1./G1b_Tb+1./Det.TES.Gt1_Tt);
Det.TES.T1o = fsolve(dE1dt,T1guess);

% Fin connector heat capacity
Det.TES.C1   =  (Det.TES.fCsn*Det.TES.gC_v.* Det.TES.T1o) * Det.TES.veff_FinCon * Det.TES.volFinCon; %J/K

% Fin connector electron-phonon coupling to bath:
Det.TES.G1b= Det.TES.nPep .* Det.TES.K1b .* Det.TES.T1o.^(Det.TES.nPep-1); %[W/K]

%Internal Thermal Conductance: evaluated at T1o 
Det.TES.Gt1_T1 = Det.TES.nt1.*Det.TES.Kt1 .* Det.TES.T1o.^(Det.TES.nt1-1);

% TES bias Power
Det.TES.Po  =   Det.TES.Ktb.* (Det.TES.Tto.^(Det.TES.nPep)-Det.fridge.T_MC.^(Det.TES.nPep)) ...
              + Det.TES.K1b.* (Det.TES.T1o.^(Det.TES.nPep)-Det.fridge.T_MC.^(Det.TES.nPep))...
              - Det.TES.Qp; %W

% TES current at equilibrium
Det.TES.Io  = sqrt(Det.TES.Po/Det.TES.Ro);%[A] 

% bias voltage needed for this equilibrium point
Det.TES.Vbias = Det.TES.Io.*(Det.elec.Rl+Det.TES.Ro);%[V]