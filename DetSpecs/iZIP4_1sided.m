function [Det]=iZIP4_1sided(fridge,elec,MatAbsb,MatTES,MatQET)
% General Detector Object for the iZIP4 in which only 1 side has been
% patterned.
%
% Inputs:
%       1) fridge: fridge object (default fUCB)
%       2) elec:  cold electronics object (default eCDMSII)
%       3) MatAbsb: material string/object for the Absorber (default 'Ge')
%       4) MatTES: material string/object for the TES (default: 'W')
%       5) MatQET: material string/object for the QP collection fins (default: 'Al')
%
% 17/05/02: MCP:  Created from iZIP4
%       
%---------------------------------------------------------
if nargin<1 | isempty(fridge)
    %if the fridge is empty, let's use the UCB fridge specs!
    fridge = fUCB();
end

if nargin<2 | isempty(elec)
    %if the electrical object is empty, let's use the CDMS-II electronics specs!
    elec = eCDMSII(fridge);
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
Det.name='iZIP4';

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single detector
Det.nP = 4;
% Number of Separately Read Out Charge Channels on a single detector
Det.nQ = 2;

%%%%%%%%%%%%%%%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The absorber object inherits all of the properties of the element object

% Absorber Dimensions
Absb.h = 25.4e-3;%[m]
Absb.r = 1.5*25.4e-3;%[m]
Absb.vol  = pi* Absb.r.^2 .* Absb.h;  %[m^3]
Absb.SAface = (pi.*Absb.r.^2); %[m^2]
Absb.SA   =  (2*pi*Absb.r).* Absb.h + 2*Absb.SAface; %[m^2];
Absb.mass = Absb.rhoD_m .* Absb.vol;

% Let's calculate the ballistic scattering length:
Absb.lscat = ScatteringLength_Ballistic(Absb.vol,Absb.SA); %[m]

%%%%%%%%%%%%%%%% TES specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TES inherited properties from its material

%number of TES/phonon channel-----
TES.n = 455;

%TES Dimensions
TES.t= 40e-9; %[m]
TES.l= 220e-6; %[m]
TES.w = 3e-6; %[m]

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
%The Aluminum collection fins are made from Al
QET=MaterialProperties('Al');

% 1D Diffusion Geometry
QET.lgc_1Ddiff = true;

QET.Afin = 142090e-12;%(m^2) area of iZIP4 QET assuming all Al is active ... a bit of an exageration
QET.lfin = 280e-6;%(m) average fin length for iZIP4
%QET.loverlap = 3.5e-6; %[m]
QET.loverlap = 5e-6; %[m]
QET.hfin=350e-9;%[m]

%let's check:
% We add 20um on each side to the TES length to account for the fact that
% the fins are wider than the TES
QET.Afin_chk = 2* QET.lfin .* (TES.l+2*20e-6);%[m^2]

%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.fSA_QET = (Det.nP*TES.n*QET.Afin)./Absb.SA;

%----- QP collection probability ------
if QET.lgc_1Ddiff
    [QET.eQPabsb,QET.ld,QET.la]= Effqp_1D_moffatt(QET.lfin,QET.hfin,QET.loverlap);
    [QET.eQPabsb_old]=Effqp_Ideal_rfin(QET.lfin, QET.ldiffQP,false);
else
    [QET.eQPabsb,QET.ld,QET.la]= Effqp_2D_moffatt(QET.lfin,QET.hfin,QET.loverlap,TES.l);
end

%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.SAactive = Det.nP*TES.n*QET.Afin;%[m^2]
%%%%%%%%%%%%%     Passive Aluminum Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The total non QET Al on the iZIP4 surface (charge bias rails/ phonon bias
% rails/bond pads/safety pads) is:

%----- rails -----
lrail = 4.1e-1;%[m]-> we almost have a meter of rails per side on the iZIP4 surface
wrail = 2*20e-6+40e-6;%[m]-> phonon rails are 20um , charge rail is 40um
SArail = lrail*wrail; %[m]

%----- safety pads -----
nsafety =369;
Apad =17671e-12;%[m^2]

SAsafety = Apad*nsafety; %[m^2]

Det.SApassive = (SArail+SAsafety); %[m^2]

% let's caclulate the fraction surface area which is covered with phonon
% absorbing Al
Det.fSA_QPabsb = (Det.SApassive+Det.SAactive)./Absb.SA;
Det.fSA_QPabsb_face = (Det.SApassive+Det.SAactive)./(Absb.SAface);

% let's calculate the percentage Al surface area which is QET and can thus
% produce signal.
Det.ePcollect = Det.SAactive/(Det.SApassive+Det.SAactive);

%%%%%%%%%%%%%     Ballistic Phonon Absorption Time     %%%%%%%%%%%%%%%%%%%%
if strcmpi(Absb.name,'Ge')
    % This was measured on G23R
    Det.tau_Pabsb= 875e-6; %[s]
elseif strcmpi(Absb.name,'Si')
    % This was hasn't been measured directly ... 
    % let's assume it's just the value on on S12Cx2
    Det.tau_Pabsb =175e-6*2; %[s]
elseif strcmpi(Absb.name,'Al2O3')
    % This has never been measured (no iZIP4 Al2O3 device!). Let's assume
    % it scales with the debye temperature:
    
    Si= MaterialProperties('Si');
    Ge= MaterialProperties('Ge');
    
    %let's use the Si ratio
    %dwdT= mean([1/750e-6/Ge.Td,1/175e-6/Si.Td]);
    dwdT= 1/175e-6/Si.Td;
    
    w= dwdT.*Absb.Td;
    
    Det.tau_Pabsb =1/w; %[s]          
elseif strcmpi(Absb.name,'ZnWO4')
    Si= MaterialProperties('Si');
    Ge= MaterialProperties('Ge');
    
    %let's use Ge as the ratio
    %dwdT= mean([1/750e-6/Ge.Td,1/175e-6/Si.Td]);
    dwdT= 1/750e-6/Ge.Td;
    
    w= dwdT.*Absb.Td;
    
    Det.tau_Pabsb =1/w; %[s]      
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
% Det.eE156 = iZIP4.eE156;
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

% Let's sum all the inductances to find the total inductance
elec.Lt = elec.Lsquid+elec.Lp+TES.L; %[H]

%//////////////// Other Measured Quantities /////////////////////
% Baseline energy resolution
Det.ResPt =[];%[sigma eVt]

%%%%%%%%%%%%%%%% Combine Objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Det.Absb=Absb;
Det.TES=TES;
Det.QET=QET;
Det.elec=elec;
Det.fridge=fridge;

