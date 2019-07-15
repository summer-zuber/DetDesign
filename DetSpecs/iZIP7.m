function [Det]=iZIP7(fridge,elec,MatAbsb,MatTES,MatQET)
% General Detector Object for the iZIP7
%
% All #'s taken from https://confluence.slac.stanford.edu/display/CDMS/Fourinch_iZIP_MIT
%
% Inputs:
%       1) fridge: fridge object (default fSNOLAB)
%       2) elec: electronics object (default eSNOLAB)
%       3) MatAbsb: material string/object for the Absorber (default 'Ge')
%       4) MatTES: material string/object for the TES (default: 'W')
%       5) MatQET: material string/object for the QP collection fins (default: 'Al')
%
% Created MCP (12/3/13)
%---------------------------------------------------------
if nargin<1 | isempty(fridge)
    %if the fridge is empty, let's use the SNOLAB fridge specs!
    fridge = fSNOLAB();
end

if nargin<2 | isempty(elec)
   elec=eSNOLAB(fridge);
end

% Initialize Absorber Object
if nargin<3 | isempty(MatAbsb)
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
if nargin<4 | isempty(MatTES)
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
if nargin<5 | isempty(MatQET)
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
Det.name='iZIP7';

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single detector
Det.nP = 12;
% Number of Separately Read Out Charge Channels on a single detector
Det.nQ = 4;

%%%%%%%%%%%%%%%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The absorber object inherits all of the properties of the element object
% Absorber Dimensions
Absb.h = 33.3e-3;%[m]
Absb.r = 50e-3;%[m]
Absb.vol  = pi* Absb.r.^2 .* Absb.h;  %[m^3]
Absb.SAface = pi.*Absb.r.^2; %[m^2]
Absb.SA   =  (2*pi*Absb.r).* Absb.h + ...
           2*Absb.SAface; %[m^2];
Absb.mass = Absb.rhoD_m .* Absb.vol;

% Let's calculate the ballistic scattering length:
Absb.lscat = ScatteringLength_Ballistic(Absb.vol,Absb.SA); %[m]

%%%%%%%%%%%%%%%% TES specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The TES object inherits all properties from the material used (usually W)

%number of TES/phonon channel-----
TES.n = 1590;

%TES Dimensions
TES.t= 40e-9; %[m]
TES.l= 160e-6; %[m] 
TES.w = 2.5e-6; %[m]

% TES fin connectors
% since the temperature in the fin connector is lower than the temperature
% in the TES, the effective volume is smaller than the true volume
TES.fFinCon=1.7;

%TES volume
TES.vol= TES.n * TES.fFinCon * (TES.t*TES.l*TES.w); %[m^3]

% Normal Resistance
TES.Rn = TES.rho .* TES.l / (TES.n*TES.w*TES.t); %[Ohm] -> Expected to be 150mOhm

%----- TES: Operating Point -----
TES.fOp=1/3;
TES.Ro = TES.Rn .* TES.fOp; %[Ohm]

%----- Channel Inductance ------
% This is a completely half-ass guess ... I should calculate more carefully
TES.L = 25e-9; %[H]

%%%%%%%%%%%%%%%%%%%%%% QET Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The QP collection fins inherit properties from their material object

% 1D Diffusion Geometry
QET.lgc_1Ddiff = true;

QET.lfin = 99e-6;%[m]
QET.hfin = 600e-9;%[m]
QET.loverlap = 3.5e-6;%[m]
QET.wfin = 200e-6;%(m)
%this is the empty space separating fins
QET.wspace = 4*8e-6; %(m)

QET.Afin = 2*QET.lfin.*(QET.wfin-QET.wspace); %(m^2) area of iZIP7 QET assuming

%----- QP collection probability ------
% from the QP trapping length and the fin geometry we can calculate the
% absorption probability ... please note that the codes being used here
% don't take into account non-perfect absorption at the interface.
%
% Ideally, in the future we would switch them out!
if QET.lgc_1Ddiff
    % 1D diffusion geometry in the fin
    [QET.eQPabsb_old]=Effqp_Ideal_rfin(QET.lfin, QET.ldiffQP,false);
    [QET.eQPabsb,QET.ld,QET.la]= Effqp_1D_moffatt(QET.lfin,QET.hfin,QET.loverlap);
else
    % 2D diffusion geometry in the fin
    %
    % for the 2D circular diffusion design, let's assume a total W/Al
    % overlap length which is the length of the TES (1/2 W/Al coverage for every fin)
    %[QET.eQPabsb]=Effqp_absb_cfin(QET.lfin,TES.ltes,QET.ldiffQP,false);
    [QET.eQPabsb,QET.ld,QET.la]= Effqp_2D_moffatt(QET.lfin,QET.hfin,QET.loverlap,TES.l);
end    

%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.SAactive = Det.nP*TES.n*QET.Afin;%[m^2]

%%%%%%%%%%%%%     Passive Aluminum Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For right now let's just assume that the passive Al coverage is solely
% due to the charge bias rails
w=40e-6; %[m]
Det.SApassive = Det.nP * TES.n .* (w * QET.lfin); %[m^2]

%------------------ Cross Check -------------------------------------------
% Estimates for passive Al coverage are coming out much larger in the lZIP
% and ulZIP ... even though naively they should be coming out much less.
% let's do some cross checks:

% width between charge and phonon lines
wQP = 1.6e-3;%[m]

%the mask has a photolithography safety edge of 2mm
w_safety = 2e-3; %[mm]
Absb.SAmask = pi*(Absb.r-w_safety).^2;%[m]

% average length of charge/phonon rails per channel
lrail = (2*Absb.SAmask/Det.nP)/(2*wQP); %[m]

% required length of a FIN to be closed packed:
QET.lfin_chk = lrail/TES.n; %[m]

%--------------------------------------------------------------------------

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

% Let's sum all the inductances to find the total inductance
elec.Lt = elec.Lsquid + elec.Lp + TES.L; %[H]

%//////////////// Measured Quantities /////////////////////
% Baseline energy resolution
Det.ResPt =[];%[sigma eVt]

%%%%%%%%%%%%%%%% Combine Objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Det.Absb=Absb;
Det.TES=TES;
Det.QET=QET;
Det.elec=elec;
Det.fridge=fridge;
