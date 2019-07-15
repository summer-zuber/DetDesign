function elec = eSLAC_G115(fridge)
% We're running G115 at SLAC with mostly SuperCDMS SNOLAB electronics ...
% but with millmax pins etc. Furthermore, the tower is all at mixing
% chamber temperature
%
% Inputs: 1) fridge object 
%
% Created MCP 12/3/13
%-------------------------------

if nargin < 1    
    fridge=fSNOLAB();
end

elec=[];

%------ Shunt Resistor ------
elec.Rs=5e-3; %[Ohm]
%elec.Rs=1e-5; %[Ohm]
elec.Ts='T_MC';

%------ Parasitic Resistances ------
elec.Rp=6e-3; %[Ohm]  Parasitic Resistance on the mixing chamber stage
%elec.Rp=1e-5; %[Ohm] check
elec.Tp={'T_MC'};

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
elec.Lsquid =25e-9; %[H]
elec.Lp = 25e-9; %[H]

%------ Current Noise ----------
elec.Si_SQUID = 4.0e-12; %[A/rtHz]
