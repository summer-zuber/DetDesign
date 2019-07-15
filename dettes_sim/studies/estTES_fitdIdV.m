function TES = estTES_fitdIdV(fitdIdV);
%Here we estimate TES parameters from fitted parameters in a dIdV step
%function
%
% input:  fitdIdV
%           .dIdV0 =
%           .w_I
%           .wp_p
%           .wp_m
%           .Rl
%           .Ro
%
% output: TES
%           .Ro
%           .Rl
%           .beta
%           .LG
%           .Lt
%           .tau0
%

% Let's invert Irwin's equations to find the core TES parameters
TES =fitdIdV;

TES.w_I = TES.w_I;
TES.w_el = TES.wp_p+TES.wp_m-TES.w_I;

E= ((TES.w_el-TES.w_I).^2-(TES.wp_p - TES.wp_m).^2)./TES.w_el./TES.w_I;

TES.beta = 4/(E+4)/TES.dIdV0/TES.Ro-TES.Rl/TES.Ro-1;
TES.LG   = (1/TES.dIdV0- (TES.Rl+ TES.Ro.*(1+TES.beta)))./(TES.Ro-TES.Rl+1/TES.dIdV0);
TES.Lt   = (TES.Rl+ TES.Ro.*(1+TES.beta))./TES.w_el;
TES.tau0 = (1-TES.LG)/TES.w_I;

% To check the above we can go backwards

%TES.tau_etf_simp = TES.tau0 ./(1+ (1-TES.Rl./TES.Ro)./(1+ TES.beta+ TES.Rl/TES.Ro).*TES.LG);%[s]
%TES.w_etf_simp   = 1./TES.tau_etf_simp;

%wp_avg = (1./TES.tau_el)/2+(1./TES.tau_I)/2;
%dw= sqrt(((1./TES.tau_el)-(1./TES.tau_I)).^2 - 4.* (TES.Ro ./elec.Lt) .* TES.LG .*(2 + TES.beta)/TES.tau0)/2;%[1/s]

%TES.wp_p = wp_avg+dw; % High Frequency Pole ->  ~ L/R
%TES.wp_m = wp_avg-dw; % Low Frequency Pole ->  ~ tau_eff

%TES.taup_p = 1./TES.wp_p;
%TES.taup_m = 1./TES.wp_m;

%TES.dIdV0 = (1-TES.LG)./(TES.Rl+TES.Ro.*(1+TES.beta)+ TES.LG.*(TES.Ro-elec.Rl));%[1/Ohm]
