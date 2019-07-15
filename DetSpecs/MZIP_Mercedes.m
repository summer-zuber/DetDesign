function [Det]=MZIP_Mercedes(fridge,MatAbsb,MatTES,MatQET,ltes,lfin,hfin,loverlap)
% General MZIP Detector Object
%
% Inputs:
%       1) fridge: fridge object (default fSNOLAB)
%       2) MatAbsb: material string/object for the Absorber (default 'Ge')
%  
%       3) MatTES: material string/object for the TES (default: 'W')
%       4) MatQET: material string/object for the QP collection fins (default: 'Al')
%       5) ltes: default = 200um
%       6) lfin: default = 200um
%       7) hfin: dafault = 600nm
%       8) loverlap: default=20um
%               
%
% Created MCP (12/3/13)
%---------------------------------------------------------
if nargin<1 | isempty(fridge)
    %if the fridge is empty, let's use the SNOLAB fridge specs!
    %fridge = fSNOLAB();
    fridge = fUCB();
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
%let's drop the TES Tc
TES.Tc = 51e-3;%[mK]

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

% default HVZIP values:
if nargin<5 | isempty(ltes)
    ltes=250e-6;%[m]
end

if nargin<6 | isempty(lfin)
    lfin=377e-6;%[m]
end

if nargin<7 | isempty(hfin)
    %hfin=600e-9;%[m]
    hfin=300e-9;%[m]
end

if nargin<8 | isempty(loverlap)
    %loverlap=3.5e-6;%[m]
    loverlap=10e-6;%[m]
end

%Initialize Detector object
Det=[];
Det.name='MZIP';

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single detector
Det.nP = 4;

%%%%%%%%%%%%%%%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The absorber object inherits all of the properties of the material object
% + ...

% Absorber Dimensions
Absb.h = 25.4e-3;%[m]

% We're going to assume that a 4" crystal is pixelized in this manner:
%   - A 1/2" thick outer active veto
%   - The 3" Boule is split into 10 same volume pixels
Absb.r = 1.5*25.4e-3;%[m]

% TES patterning isn't possible very near the outer edge
Absb.w_safety = 2e-3; %[m]

Absb.vol  = pi* Absb.r.^2 .* Absb.h;  %[m^3]
Absb.SAface = pi*Absb.r .^2; %[m^2]

Absb.SApattern = pi*(Absb.r-Absb.w_safety).^2;%[m^2]
Absb.SA   =  2*Absb.SAface+ (2*pi*Absb.r).*Absb.h;%[m^2]

% Absorber Mass
Absb.mass = Absb.rhoD_m .* Absb.vol;

% Let's calculate the ballistic scattering length:
Absb.lscat = ScatteringLength_Ballistic(Absb.vol,Absb.SA); %[m]

%%%%%%%%%%%%%%%% TES specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TES.t = 40e-9; %[m] TES thickness
TES.l = ltes;%[m]   TES length
TES.w = 2.5e-6; %[m]TES width

% TES fin connectors
% since the temperature in the fin connector is lower than the temperature
% in the TES, the effective volume is smaller than the true volume
TES.fFinCon=3.0;

%---------------- Calculate Resistance ------------------------------------
% Due to electronics constraints, the normal resistance should be ~300mOhm

% the resistance of a single TES is
Rn_1tes = TES.rho .* TES.l / (TES.w*TES.t); %[Ohm]

if false
    %Rn_wanted = 200e-3; %[Ohm]
    Rn_wanted = 400e-3; %[Ohm]

    % Therefore, the number of TES which can be put in parallel is
    TES.n = ceil(Rn_1tes/Rn_wanted);

    % Now let's recalculate the Rn value, just a bit more precisely
    TES.Rn = Rn_1tes/TES.n; %[Ohm]
else
    TES.n = 889;    
    TES.Rn = Rn_1tes/TES.n; %[Ohm]
end    

%--------------- W Volume -------------------------------------------------
%TES volume
TES.vol= TES.n * TES.fFinCon * (TES.t*TES.l*TES.w); %[m^3]

%----- TES: Operating Point -----
TES.fOp = 1/3;
TES.Ro  = TES.Rn .* TES.fOp; %[Ohm]

%----- Channel Inductance ------
% This is a completely half-ass guess ... I should calculate more carefully
TES.L = 25e-9; %[H]

%%%%%%%%%%%%%%%%%%%%%% QET Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The collection fin object inherits information from its material object
% 2D Diffusion Geometry
QET.lgc_1Ddiff = false;

QET.lfin = lfin;%[m]
QET.hfin = hfin; %[m]
QET.loverlap = loverlap; %[m]

if QET.lgc_1Ddiff
    [QET.eQPabsb,QET.ld,QET.la]= Effqp_1D_moffatt(QET.lfin,QET.hfin,QET.loverlap);
else
    [QET.eQPabsb,QET.ld,QET.la]= Effqp_2D_moffatt(QET.lfin,QET.hfin,QET.loverlap,TES.l);
end

%------ QET Active Area ----- 

% number of fins in a QET
if QET.lfin > 100e-6
    QET.nfin = 12; %# of fins in a QET
else
    QET.nfin = 4;
end

if false
    %wempty = 6e-6;%[m]
    wempty = 8e-6;%[m]
    QET.Afin_empty = (QET.nfin*QET.lfin+TES.l)*wempty; %[m^2]
    QET.Afin = pi* QET.lfin .*(QET.lfin + TES.l/2) - QET.Afin_empty; %[m^2]
    %QET.Afin = pi* QET.lfin.^2 + 2*QET.lfin*TES.l - QET.Afin_empty; %[m^2]
else
    % Pulled directly from mask design file
    QET.Afin= 521940e-12;%[m^2]
end
    
%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.SAactive = Det.nP*TES.n*QET.Afin;%[m^2]

%%%%%%%%%%%%% Passive Aluminum Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pulled Directly from mask TDB file
Det.SApassive = 353361378e-12; %[m^2]

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
               .* DiZIP4.fSA_QPabsb./Det.fSA_QPabsb;
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
