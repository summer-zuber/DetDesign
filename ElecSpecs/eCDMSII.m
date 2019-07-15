function elec = eCDMSII(fridge)
% Creates an object holding all of the specifications for the old CDMS-II electronics
%
% Inputs: 1) fridge object (Default: Soudan Fridge)
%
% Created MCP 12/3/13
%-------------------------------

if nargin < 1    
    fridge=fSoudan();
end

elec=[];

%------ Shunt Resistor ------
elec.Rs=24e-3; %[Ohm]
elec.Ts='T_Still';

%------ Parasitic Resistances ------
elec.Rp=[6e-3  %[Ohm] Parasitic Resistance on the still stage
         2e-3  %[Ohm] Parasitic Resistance on the 4K stage
         4e-3];%[Ohm] Parasitic Resistance on the mixing chamber stage

% elec.Rp=[4e-3  %[Ohm] Parasitic Resistance on the still stage
%          2e-3  %[Ohm] Parasitic Resistance on the 4K stage
%          2e-3];%[Ohm] Parasitic Resistance on the mixing chamber stage
     
elec.Tp={'T_Still'
         'T_4K'
         'T_MC'};

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
elec.Lsquid =250e-9; %[H]
elec.Lp = 100e-9; %[H]

%------ Current Noise ----------
elec.Si_SQUID = 2e-12; %[A/rtHz]
%elec.Si_SQUID= 8e-12;%[A/rtHz]
