function [Det]=ulZIP(fridge,MatAbsb,MatTES,MatQET)
% General Detector Object for the ultra light mass WIMP detector
%
% Inputs:
%       1) fridge: fridge object (default fSNOLAB)
%       2) MatAbsb: material string/object for the Absorber (default 'Ge')
%       3) MatTES: material string/object for the TES (default: 'W')
%       4) MatQET: material string/object for the QP collection fins (default: 'Al')
%
% Created MCP (12/3/13)
%---------------------------------------------------------
if nargin<1 | isempty(fridge)
    %if the fridge is empty, let's use the SNOLAB fridge specs!
    fridge = fSNOLAB();
end

% Initialize Absorber Object
if nargin<2 | isempty(MatAbsb)
    % if there is no input value, let's assume that the absorber is Ge
    Absb = MaterialProperties('Ge');
elseif iscellstr(MatAbsb)| ischar(MatAbsb)
    Absb = MaterialProperties(MatAbsb);
else
    % let's inherit all the physical properties from the absorber material
    % to our absorber object
    Absb = MatAbsb;
end  

%Initialize TES Object
if nargin<3 | isempty(MatTES)
    % if there is no input value, let's assume that the TES is W
    TES = MaterialProperties('W');
elseif iscellstr(MatTES)| ischar(MatTES)
    TES = MaterialProperties(MatTES);
else
    % let's inherit all the physical properties from the TES material
    % to our TES object
    TES = MatTES;
end  

%Initialize QET Object
if nargin<4 | isempty(MatQET)
    % if there is no input value, let's assume that the QET fins are made
    % of Al
    QET = MaterialProperties('Al');
elseif iscellstr(MatQET)| ischar(MatQET)
    QET = MaterialProperties(MatQET);
else
    % let's inherit all the physical properties from fin material
    % to our absorber object
    QET = MatQET;
end  

%Initialize Detector object
Det=[];
Det.name='ulZIP';

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single pixel
Det.nP = 1;
% Number of Separately Read Out Charge Channels on a single detector
Det.nQ = 0;

%%%%%%%%%%%%%%%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The absorber object inherits all of the properties of the material object
% + ...

% Absorber Dimensions
Absb.h = 33.3e-3;%[m]

% We're going to assume that a 4" crystal is pixelized in this manner:
%   - A 1/2" thick outer active veto
%   - The 3" Boule is split into 10 same volume pixels
Absb.r = 37.5e-3;%[m]
Absb.vol  = pi* Absb.r.^2 .* Absb.h /10;  %[m^3]

% The precise surface area depends upon the pixelization geometry 
% (and probably will be slightly different for the various pixels)
% Let's assume rectangular crystal geometries
d = sqrt(Absb.vol/Absb.h); %[m]
Absb.SA   =  2*d.^2+4*d*Absb.h;  %[m^2];

% Absorber Mass
Absb.mass = Absb.rhoD_m .* Absb.vol;

% Let's calculate the ballistic scattering length:
Absb.lscat = ScatteringLength_Ballistic(Absb.vol,Absb.SA); %[m]

%%%%%%%%%%%%%%%% TES specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TES properties inherited from material object

% For this device we want the smallest possible TES volume that still has
% reproducible transitions. One of the problems with going to shorter and
% shorter TES is that eventually proximity effect will increase the Tc too
% much.

% So let's use 40um as our lower edge:
TES.t=40e-9; %[m]

%let's have 2 different designs
if MatQET.ldiffQP < 150e-6
    TES.l= 80e-6; %[m] % 1.5ms for 135e-6 QP trapping length
else
    TES.l = 25e-6; %[m] % 1.3ms for 300e-6 QP trapping length
end

TES.w=2.5e-6; %[m]

% TES fin connectors
% since the temperature in the fin connector is lower than the temperature
% in the TES, the effective volume is smaller than the true volume
TES.fFinCon=1.7;

%---------------- Calculate Resistance ------------------------------------
% Due to electronics constraints, the normal resistance should be ~150mOhm
Rn_wanted = 150e-3; %[Ohm]

% while the resistance of a single TES is
Rn_1tes = TES.rho .* TES.l / (TES.w*TES.t); %[Ohm]

% Therefore, the number of TES which can be put in parallel is
TES.n = ceil(Rn_1tes/Rn_wanted);

% Now let's recalculate the Rn value, just a bit more precisely
TES.Rn = Rn_1tes/TES.n; %[Ohm]

%--------------- W Volume -------------------------------------------------
%TES volume
TES.vol= TES.n * TES.fFinCon * (TES.t*TES.l*TES.w); %[m^3]

%----- TES: Operating Point -----
TES.fOp=1/3;
TES.Ro = TES.Rn .* TES.fOp; %[Ohm]

%----- Channel Inductance ------
% This is a completely half-ass guess ... I should calculate more carefully
TES.L = 25e-9; %[H]

%%%%%%%%%%%%%%%%%%%%%% QET Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The collection fin object inherits information from its material object
% 2D Diffusion Geometry
QET.lgc_1Ddiff = false;

% Maximum Energy S/N occurs with a fin length ~75% of the trapping
% length for 2D Diffusion Geometry
QET.lfin = 0.75 .* QET.ldiffQP; %[m]

% Empty Separation Space between fins so that there are no shorts:
wempty = 8e-6;%[m]
Aempty = 4*QET.lfin*wempty; %[m^2]

%QET.Afin = pi* QET.lfin .*(QET.lfin + TES.l/2) - Aempty; %[m^2]
QET.Afin = pi* QET.lfin.^2 + 2*QET.lfin*TES.l - Aempty; %[m^2]


%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.fSA_QET = (Det.nP*TES.n*QET.Afin)./Absb.SA;

%----- QP collection probability ------
% from the QP trapping length and the fin geometry we can calculate the
% absorption probability ... please note that the codes being used here
% don't take into account non-perfect absorption at the interface.
%
% Ideally, in the future we would switch them out!
if QET.lgc_1Ddiff
    % 1D diffusion geometry in the fin
    [QET.eQPabsb]=Effqp_Ideal_rfin(QET.lfin, QET.ldiffQP,false);
else
    % 2D diffusion geometry in the fin
    %
    % for the 2D circular diffusion design, let's assume a total W/Al
    % overlap length which is x2 the length of the TES
    [QET.eQPabsb]=Effqp_Ideal_cfin(QET.lfin,2*TES.l,QET.ldiffQP,false);
end    

%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.SAactive = Det.nP*TES.n*QET.Afin;%[m^2]

%%%%%%%%%%%%% Passive Aluminum Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Since the QETs are relatively sparse, we need bias rails to travel
% between the QETs

% rail width
wrail=10e-6;%[m]

% to calculate rail length, let's take the surface area of an instrumented
% face and divide by the number of QETs to get the lot area per QET
SAlot= Absb.SA/6/TES.n; %[m^2]

% if we think of these lots as squares, then the length of the rails is 
llot=sqrt(SAlot); %[m]

% let's just multiply this by the number of QETs, to get the total length
% of rails, and then multiply by rail width to get the total rail area

Det.SApassive = TES.n * llot * wrail; %[m^2]

% let's caclulate the fraction surface area which is covered with phonon
% absorbing Al
Det.fSA_QPabsb = (Det.SApassive+Det.SAactive)./Absb.SA;

% let's calculate the percentage Al surface area which is QET and can thus
% produce signal.
Det.ePcollect = Det.SAactive/(Det.SApassive+Det.SAactive);

%%%%%%%%%%%%%     Ballistic Phonon Absorption Time     %%%%%%%%%%%%%%%%%%%%
% Since this device has not been tested/created, we need to estimate the
% ballistic phonon collection time from that of the iZIP4.

%let's start with an iZIP4 design which uses the correct crystal
DiZIP4 = iZIP4(fridge,elec,Absb.name);

% The characteristic phonon collection time is proportional to:
% 1) the absorber scattering length
% 2) the fractional surface area coverage of the QETs (inversely)
Det.tau_Pabsb = DiZIP4.tau_Pabsb ...
               .* Absb.lscat./DiZIP4.Absb.lscat ...
               .* DiZIP4.fSA_QET./Det.fSA_QET;
Det.w_Pabsb=1/Det.tau_Pabsb; %[rad/s]           
           
% if strcmpi(Absb.name,'Ge')
%     % This was measured on G48
%     Det.tau_Pabsb= 750e-6; %[s]
% elseif strcmpi(Absb.name,'Si')
%     % This was measured on S12C
%     Det.tau_Pabsb =175e-6; %[s]
% else
%     display('iZIP4 Phonon Collection Specs only available for Ge/Si')
%     return
% end 
% Det.w_Pabsb=1/Det.tau_Pabsb; %[rad/s]

%%%%%%%%%%%%%%% Total Phonon Energy Collection Efficiency %%%%%%%%%%%%%%%%%
 
% The loss mechanisms in our detector are:
% 1) subgap downconversion of athermal phonons in the crystal
% 2) collection of athermal phonons by passive metal on the surface of our
% detector ( Det.ePcollect)
% 3) Efficiency of QP production in Al fin (QET.ePQP)
% 4) Efficiency of QP transport to TES (QET.eQPabsb)
% 5) Energy conversion efficiency at W/Al interface
% 6) ?

% Let's combine 1), 5), and 6) together and assume that it is the same as the measured/derived value from iZIP4
Det.eE156 = DiZIP4.eE156;

% Then we can estimate the total phonon energy collection efficiency

%               1+5+6       2               3               4
Det.eEabsb = Det.eE156 .* Det.ePcollect .* QET.ePQP .* QET.eQPabsb;

%Rather than estimate Det.eEAbsb, let's take the measured Det.eEabsb, and
%back calculate to find Det.eE156
%Det.eEabsb = 0.13;
%Det.eE156 = Det.eEabsb / (Det.ePcollect .* QET.ePQP .* QET.eQPabsb);

%%%%%%%%%%%%%%%%%%% Thermal Conductance to Bath %%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Thermal Conductance Coefficient between Detector & Bath -----
Det.Kpb = 1.55e-4; %[W/K^4]
%----- Thermal Conductance Coefficient between Detector & Bath
Det.nKpb = 4;  

%%%%%%%%%%%%%%%% Electronics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This device uses low impedance SNOLAB electronics
elec= eSNOLAB(fridge);
% Let's sum all the inductances to find the total inductance
elec.Lt = elec.Lsquid+elec.Lp+TES.L; %[H]

%//////////////// Measured Quantities /////////////////////
% Baseline energy resolution
Det.ResPt =[];%[sigma eVt]

%%%%%%%%%%%%%%%% Combine Objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Det.Absb=Absb;
Det.TES=TES;
Det.QET=QET;
Det.elec=elec;
Det.fridge=fridge;
