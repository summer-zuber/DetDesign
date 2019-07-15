function elec = eOPD(fridge)
% Creates an object holding all of the specifications needed for the
% ultimate optical phonon detector
%
% Inputs: 1) fridge object (Default: OPD fridge)
%
% Created MCP 18/05/12
%-------------------------------

if nargin < 1    
    fridge=fOPD();
end

elec=[];

%------ Shunt Resistor ------
%elec.Rs=5e-3; %[Ohm]
elec.Rs=1e-3; %[Ohm]
elec.Ts='T_CP';

%------ Parasitic Resistances ------
%elec.Rp = [2e-3,2e-3]; %[Ohm]
elec.Rp = [.1e-3,.1e-3]; %[Ohm]
elec.Tp={'T_MC','T_CP'};

nRp= length(elec.Rp);
% let's make overall shunt quantities         
elec.Rpt = sum(elec.Rp);

elec.Tpt = 0;
for jj=1:nRp
    elec.Tpt = elec.Tpt + elec.Rp(jj).*fridge.(elec.Tp{jj});
end
elec.Tpt = elec.Tpt/elec.Rpt; %[K]

%------ Total Load Resistance ------
elec.Rl = elec.Rpt+elec.Rs; %[Ohm]
elec.Tl = (elec.Rpt*elec.Tpt + elec.Rs*fridge.(elec.Ts))/elec.Rl; %[K]

%------ Inductance -----------------
elec.Lsquid =400e-9; %[H]
elec.Lp = 25e-9; %[H]

%------ Current Noise ----------
elec.Si_SQUID = 2e-12; %[A/rtHz]

