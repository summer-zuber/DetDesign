function dydt = dydt_ParallelTES_eq(t,y,Det)
% this function computes the time derivative of the state vector for a 
% system of TES in parallel readout by the same SQUID (i.e. they share 
% Rl and Lt)
%
% inputs:
%   1) t
%   2) y= [I;T] -> where I and T are both column vectors
%   3) Det -> the detector object
%
% 3/31/17 MCP created ->  based on previous simulation code
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pc= PhysicalConstants;

%------------- Steady State Constraints ----------------------------------
%the inhomogenous terms of the NDE (heat and voltage fluctuations) are 
%assumed to be zero
dVbias   = 0;
dQ = 0;

%-------------------------------------------------------------------------
%FIRST let's split y into it's more natural variables
I = y(1:Det.TES.nSim); %[A]
T = y(Det.TES.nSim+1:end); %[K]

%Total current flowing through all the parallel TES.
Itot = sum(I); %[A]

%----- Calculate TES Resistance ---------------------------
[R,alpha,beta]= Rtes(I,T,Det);

Rparallel = Rtes_parallel(R); %[Ohm]

%----- Calculate Heat Capacity ----------------------------
C =  (Det.TES.fCsn*Det.TES.gC_v.*T)*(Det.TES.vol); %J/K

%---- Calculate Cooling Power -----------------------------
Pcool = Det.TES.gPep_v*(Det.TES.vol).* (T.^(Det.TES.nPep)-Det.fridge.T_MC.^(Det.TES.nPep)); %[W]

%----- Calculate dTdt -------------------------------------
% C dTdt = I^2 R - Pcool + dQ
dTdt = ( I.^2 .*R - Pcool + dQ)./C; %[K/s]

%----- Calculate dItotdt -----------------------------------
% Convert from current biasing to voltage biasing
Vbias = Det.TES.Ibias .* Det.elec.Rs;

dItotdt =(Vbias+dVbias-Itot.*(Det.elec.Rl+Rparallel))/Det.elec.Lt;

%----- Calculate the individual dIdt -----------------------
%let's follow the calculation done in "Parallel_TES_Dynamics.pdf"

dRparalleldT = Rparallel.^2 .* alpha./(R.*T);
dRparalleldI = Rparallel.^2 .*  beta./(R.*I);

Zj =  dItotdt.*Rparallel + ...
      Itot.*sum(dRparalleldT.*dTdt) +...
     -I.*alpha.*R.*dTdt./T;


Mjk= diag(R.*(1+beta)) - repmat((Itot.*dRparalleldI)',Det.TES.ntessim,1);

dIdt = Mjk\Zj;

%Finally let's recombine the time derivatives into  a single vector dydt
dydt=[dIdt;dTdt];