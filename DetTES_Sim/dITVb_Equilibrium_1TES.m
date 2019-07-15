function dY= dITVb_Equilibrium_1TES(X,Det)
% The non-linear equations which define the equilibrium
% point for a simple TES.
%
% Let's use these equations to find the equilibrium point and the wanted
% voltage bias:
%
% 1)  L dIdt = Vb - I (Rl +R)
% 2)  C dTdt = I^2 R - Pcool
% 3)     0   = R - Rwanted
%
% X = Equilibrium Point Phase Space Vector
%     I ->  TES current 
%     T ->  TES temperature
%     Vb ->  Bias Voltage
%
% Det:  Detector Object which contains TES properties (Tc etc.)
%
% 3/31/17 MCP created ->  based on previous simulation code
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pc= PhysicalConstants;

%-------------------------------------------------------------------------
%FIRST let's split y into it's more natural variables
I  = X(1);
T  = X(2);
Vb = X(3);

% Calculate TES Resistance

%let's calculate the time rate of change of both the TES current, I,
% and the TES temperature,T. To do this let's first find the resistance 
%for these variables Rtes =R_tes(I,T);
[R,alpha,beta]= Rtes(I,T,Det);

% We want the TES to have a resistance Ro
dR = R-Det.TES.Ro;

%Since we have R we can calculate dIdt
% use 1) voltage equation
% L dI/dt=Vbias - I*(Rl+Rtes)
dIdt = (Vb-I*(Det.elec.Rs+R))/Det.elec.Lt;

% C dTdt = I^2 Rtes - Pcool +dQtes

% Heat Capacity
C =  (Det.TES.fCsn*Det.TES.gC_v.* T)* Det.TES.vol; %J/K

dTdt = ( I^2*R + ...  % Joule Heating %[W]
        - Det.TES.gPep_v*(Det.TES.vol).* (T.^(Det.TES.nPep)-Det.fridge.T_MC.^(Det.TES.nPep)) ...%Cooling to the substrate %[W]      
       )./C; %[K/s]
   
%Let's combine all 3 equations into an output vector:
dY=[dIdt;dTdt;dR];