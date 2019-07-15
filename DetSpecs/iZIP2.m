function [Det]=iZIP2(fridge,MatAbsb,MatTES,MatQET)
% General Detector Object for the iZIP2
%
% Inputs:
%       1) fridge: fridge object (default fUCB)
%       2) MatAbsb: material string/object for the Absorber (default 'Ge')
%       3) MatTES: material string/object for the TES (default: 'W')
%       4) MatQET: material string/object for the QP collection fins (default: 'Al')
%
% Created MCP (12/3/13)
%---------------------------------------------------------
if nargin<1 | isempty(fridge)
    %if the fridge is empty, let's use the UCB fridge specs!
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
Det.name='iZIP2';

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single detector
Det.nP = 5;
% Number of Separately Read Out Charge Channels on a single detector
Det.nQ = 4;

%%%%%%%%%%%%%%%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The absorber object inherits all of the properties of the element object

% Absorber Dimensions
Absb.h = 25.4e-3;%[m]
Absb.r = 1.5*25.4e-3;%[m]
Absb.vol  = pi* Absb.r.^2 .* Absb.h;  %[m^3]
Absb.SA   =  (2*pi*Absb.r).* Absb.h + ...
           2*(pi.*Absb.r.^2); %[m^2];
Absb.mass = Absb.rhoD_m .* Absb.vol;

% Let's calculate the ballistic scattering length:
Absb.lscat = ScatteringLength_Ballistic(Absb.vol,Absb.SA); %[m]

%%%%%%%%%%%%%%%% TES specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%number of TES/phonon channel-----
TES.n = 556;

%TES Dimensions
TES.t= 40e-9; %[m]
TES.l= 220e-6; %[m]
TES.w = 2.5e-6; %[m]

% TES fin connectors
% since the temperature in the fin connector is lower than the temperature
% in the TES, the effective volume is smaller than the true volume
TES.fFinCon=1.7;

%TES volume
TES.vol= TES.n * TES.fFinCon * (TES.t*TES.l*TES.w); %[m^3]

% Normal Resistance
TES.Rn = TES.rho .* TES.l / (TES.n*TES.w*TES.t); %[Ohm]

%----- TES: Operating Point -----
TES.fOp=1/3;
TES.Ro = TES.Rn .* TES.fOp; %[Ohm]

%----- Channel Inductance ------
% This is a completely half-ass guess ... I should calculate more carefully
TES.L = 100e-9; %[H]

%%%%%%%%%%%%%%%%%%%%%% QET Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The QET Fin object inherited properties from it's material

% 1D Diffusion Geometry
QET.lgc_1Ddiff = true;

QET.Afin = 158538e-12;%(m^2) area of iZIP4 QET assuming all Al is active ... a bit of an exageration
QET.lfin = 300e-6;%(m) average fin length for iZIP4

%let's check:
% We add 20um on each side to the TES length to account for the fact that
% the fins are wider than the TES
QET.Afin_chk = 2* QET.lfin .* (TES.l+2*20e-6);%[m^2]

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
    % overlap length which is the length of the TES (1/2 W/Al coverage for every fin)
    [QET.eQPabsb]=Effqp_absb_cfin(QET.lfin,TES.ltes,QET.ldiffQP,false);
end    

%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.SAactive = Det.nP*TES.n*QET.Afin;%[m^2]

%%%%%%%%%%%%%     Passive Aluminum Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%let's just assume that the passive Al is %0.5 of the total surface area

Det.SApassive = 0.005*Absb.SA; %[m^2]

% let's caclulate the fraction surface area which is covered with phonon
% absorbing Al
Det.fSA_QPabsb = (Det.SApassive+Det.SAactive)./Absb.SA;

% let's calculate the percentage Al surface area which is QET and can thus
% produce signal.
Det.ePcollect = Det.SAactive/(Det.SApassive+Det.SAactive);

%%%%%%%%%%%%%     Ballistic Phonon Absorption Time     %%%%%%%%%%%%%%%%%%%%
if strcmpi(Absb.name,'Ge')
    % This was measured on G3D
    Det.tau_Pabsb= 850e-6; %[s]
elseif strcmpi(Absb.name,'Si')
    % We're going to assume an identical Si/Ge scaling factor as seen in
    % iZIP4 detectors
    Det.tau_Pabsb =(175/750) * 850e-6; %[s]
else
    display('iZIP4 Phonon Collection Specs only available for Ge/Si')
    return
end
Det.w_Pabsb=1/Det.tau_Pabsb; %[rad/s]

%%%%%%%%%%%%%%% Total Phonon Energy Collection Efficiency %%%%%%%%%%%%%%%%%
 
% The loss mechanisms in our detector are:
% 1) subgap downconversion of athermal phonons in the crystal
% 2) collection of athermal phonons by passive metal on the surface of our
% detector ( Det.ePcollect)
% 3) Efficiency of QP production in Al fin (QET.ePQP)
% 4) Efficiency of QP transport to TES (QET.eQPabsb)
% 5) Energy conversion efficiency at W/Al interface
% 6) ?

% % Let's combine 1), 5), and 6) together and assume that it is the same as the measured/derived value from iZIP4
% Det.eE156 = DiZIP4.eE156;
% 
% % Then we can estimate the total phonon energy collection efficiency
% 
% %               1+5+6       2               3               4
% Det.eEabsb = Det.eE156 .* Det.ePcollect .* QET.ePQP .* QET.eQPabsb;

%Rather than estimate Det.eEAbsb, let's take the measured Det.eEabsb, and
%back calculate to find Det.eE156
Det.eEabsb = 0.13;
Det.eE156 = Det.eEabsb / (Det.ePcollect .* QET.ePQP .* QET.eQPabsb);

%%%%%%%%%%%%%%%%%%% Thermal Conductance to Bath %%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Thermal Conductance Coefficient between Detector & Bath -----
Det.Kpb = 1.55e-4; %[W/K^4]
%----- Thermal Conductance Coefficient between Detector & Bath
Det.nKpb = 4;  

%%%%%%%%%%%%%%%% Electronics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This device uses CDMS-II electronics
elec= eCDMSII(fridge);
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

