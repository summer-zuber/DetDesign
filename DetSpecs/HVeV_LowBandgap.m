%Low Bandgap Semiconductors

InSb       = [];
InSb.Eg    =;
InSb.n_n   = 5e13% [dopants/cm^3]
InSb.rho_m = ;
InSb.Tdebye= 203;%[K]

%supplier: https://www.mtixtl.com/PbS-a-101005S1-3.aspx
PbSe            = [];
PbSe.Eg         = 0.27; %[eV]
PbSe.rho_m      = 8.15;%[g/cm^3]
PbSe.harndness  = 2.7;%[mohs]
PbSe.Tdebye     = 145;%[145-225K]

%supplier: https://www.mtixtl.com/PbS-a-101005S1-3.aspx
%https://www.powerwaywafer.com/pbte-single-crystal-substrate.html
% Seems to have an enormous amount of anharmonic decay :(
PbTe       = [];
PbTe.Eg    = 0.32;%[eV]
PbTe.rho_m = 8.16;%[g/cm^3]
PbTe.harndness = 3;%[mohs]


