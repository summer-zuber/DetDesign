function [Det] = SimpleEquilibrium_1TES(Det,beta,Qp)
% A Quick/Inaccurate Estimate the TES equilibrium point. Basically, we
% calculate the bias point assuming that beta is equal to zero, and then we
% adhoc add a beta after the fact.
%
% It's a bit sketchy since it doesn't fit the functional form that we assume for TES resistance,
% However, it will produce an dynamical system that will exactly match all small signal measurements of a simple 1 block TES
%   
%  Assumptions:
%       1) beta = 0;f
%       2) R(T)= Rn(1+tan((T-Tc)/w_tc))/2

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

%let's check to see if the Det object is carrying a field for DC parasitic
%power loading


if nargin <3 || isempty(Qp)  
    %If the detector object is carrying a parasitic heating variable, then
    %we just leave it ... if not let's add a field and set it to zero
    if ~isfield(Det.TES,'Qp')
        Det.TES.Qp=0;
    end 
else
    %let's set the 
    Det.TES.Qp=Qp;%[W]
end

%\\\\\\\\\\\ Physical Constants and Invariant Characteristics \\\\\\\\\\\\\
pc=PhysicalConstants();

%------- Bias Point Temperature ------------------------------------
% let's calculate the temperature of the operating point resistance.
% [Notice that if the resistance of the TES changes with current, then this doesn't work]

zeta_o= log(Det.TES.fOp/(1-Det.TES.fOp))/2;
Det.TES.To = zeta_o * Det.TES.wTc + Det.TES.Tc; %K


%------- Alpha/Beta at Transition Point -----------------------------------

% let's calculate bias location along the normalized transition curve

%let's also calculate alpha 
Det.TES.alpha= 2.*Det.TES.To/Det.TES.wTc ./ exp(zeta_o) ./ (exp(zeta_o)+exp(-zeta_o));
Det.TES.beta = beta;

%------ TES Properties @ Equilibrium --------------------------------------
% TES heat capacity
Det.TES.C   =  (Det.TES.fCsn*Det.TES.gC_v.* Det.TES.To)*(Det.TES.vol); %J/K

% TES phonon/electron thermal coupling 
Det.TES.Gep =  Det.TES.nPep*Det.TES.gPep_v*(Det.TES.vol).*  Det.TES.To.^(Det.TES.nPep-1); %W/K

% TES bias Power
Det.TES.Po  =               Det.TES.gPep_v*(Det.TES.vol).* (Det.TES.To.^(Det.TES.nPep)-Det.fridge.T_MC.^(Det.TES.nPep))-Det.TES.Qp; %W

%loop Gain
Det.TES.LG  = Det.TES.alpha * Det.TES.Po ./ (Det.TES.Gep .* Det.TES.To);
%if Det.TES.LG<0
%    Det.TES.LG=0;
%end    

% TES current at equilibrium
Det.TES.Io  = sqrt(Det.TES.Po/Det.TES.Ro);%[A] 

% bias voltage needed for this equilibrium point
Det.TES.Vbias = Det.TES.Io.*(Det.elec.Rl+Det.TES.Ro);%[V]

% simple thermal conductance
Det.TES.tau0= Det.TES.C./Det.TES.Gep; %[s]

%------ Sensor Bandwidth --------------------------------------------------
%Let's calculate the sensor bandwidth
Det.TES.tau_etf= Det.TES.tau0./(1+Det.TES.LG.*(1-Det.elec.Rl./Det.TES.Ro)./(1+Det.TES.beta+Det.elec.Rl./Det.TES.Ro));
Det.TES.w_etf= 1./Det.TES.tau_etf;
