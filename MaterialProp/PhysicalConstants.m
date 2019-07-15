function pc= PhysicalConstants();
%this is an object containing various physical objects (most in SI units)

pc=[];

%boltzman constant
pc.kb    = 1.3806503e-23; %[J/K]

%Stefan-Boltzman Black Body Constant
pc.sigma_SB= 5.670367e-8;%[W/m^2/K^4]

%electronic charge in vacuum
pc.qe    = 1.60219e-19; % [C]

%Avagardo's number
pc.nA_g  = 6.02214129e23; %[Atoms/gmol] -> Avagardo's number in g
pc.nA_kg = pc.nA_g*1e3;% [Atoms/kmol] -> Avagardo's number in kg

%hbar
pc.hbar= 1.05459e-34;%[Js]
pc.h = pc.hbar*(2*pi);%[Js]

%speed of light
pc.c=3e8; %m/s

%mass of electron
pc.Me_kg=511*1e3*pc.qe/pc.c.^2; %kg

% Lorenz number
% proportionality coefficient in the weiderman-franz law
% K= Lwf T/rho
pc.Lwf=2.44e-8; %[W Ohm / K^2]

%vacuum permitivity
pc.epsilon0 = 8.854187e-12; %[F/m]

%Gravitational Constant
pc.G = 6.674e-11;%[m^3/(kg s^2)]