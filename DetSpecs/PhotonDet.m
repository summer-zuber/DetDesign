function [Det]=PhotonDet(fridge,MatAbsb,MatTES,MatQET,ltes,lfin,hfin,loverlap)
% General Detector Object for the light mass WIMP detector (usually with luke neganov gain)
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
    %fridge = fUCB(); 
    %fridge=fOPD();
    
end
%This device uses low impedance SNOLAB electronics
%elec= eSNOLAB(fridge);
%elec= eCDMSII(fridge);
elec= eSLAC_G115(fridge);
%elec= eOPD(fridge);

% Initialize Absorber Object
if nargin<2 | isempty(MatAbsb)
    % if there is no input value, let's assume that the absorber is Ge
    Absb = MaterialProperties('Si');
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
%TES.Tc=41.5e-3; %[K]

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

% default Photon Detector values:
if nargin<5 | isempty(ltes)
    ltes=140e-6;%[m]
end

if nargin<6 | isempty(lfin)
    lfin=200e-6;%[m]
end

if nargin<7 | isempty(hfin)
    hfin=600e-9;%[m]
end

if nargin<8 | isempty(loverlap)
    % we designed 2 W TES masks:
    %1) 10um W/Al overlap used on wafer #10 (PD1) and #13 (PD2
    loverlap=10e-6;%[m]  
    
    %2) 20um W/Al overlap used on wafer #9
    %loverlap=20e-6;%[m] 
end

lgc_build = true;

%Initialize Detector object
Det=[];
Det.name='Photon Det';

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single detector
Det.nP = 1;

%%%%%%%%%%%%%%%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The absorber object inherits all of the properties of the material object
% + ...

% Absorber Dimensions
Absb.h = 1e-3;%[m]
Absb.r = 38.1e-3;%[m]
% TES patterning isn't possible very near the outer edge
Absb.w_safety = 3e-3; %[m]

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

if lgc_build
    TES.w = 4e-6; %[m]TES width
else
    TES.w = 1.5e-6; %[m]TES width
end

% number of fins in a QET
if lfin > 100e-6
    QET.nfin = 6; %# of fins in a QET
else
    QET.nfin = 4;
end    

% fraction of the Al fin edge that is adjacent to the TES which is covered with Al 
TES.foverlap_width= 1.2;

% Volume of a single TES
TES.vol1_TES= TES.w * TES.l * TES.t; %[m^3]

% Volume of the W only Fin connector
% orginal guess
TES.vol1_WFinCon = QET.nfin*2.5e-6*4e-6*TES.t ... % [m^3]the neck of the W only fin connector
                   +2.5e-6*(2*TES.l*TES.foverlap_width)*TES.t;% [m^3]the long piece W only piece going along the fin

% Looked in L-EDIT for exact dimensions 
TES.area1_WFinCon_midchk =(6e-6+24e-6)*2.5e-6+3e-6*4.4e-6+1e-6*1e-6;%[m^2]
TES.area1_WFinCon_mid = 88.905e-12;%[m^2]

TES.area1_WFinCon_endchk=7.5e-12+22.25e-12+37.5e-12+1e-12;%[m^2]
TES.area1_WFinCon_end=68.390e-12;%[m^2]

TES.vol1_WFinCon_LEDIT= (2*TES.area1_WFinCon_end+ 4*TES.area1_WFinCon_mid)*TES.t; %[m^3]
 
% Volume of the W/Al overlap                            
TES.vol1_WAloverlap = loverlap*(2*TES.l*TES.foverlap_width)*TES.t;

%Looked in L-EDIT For exact dimensions:
if true
   % wafer #13(PD2) and wafer #10 (PD1) have 10um overlap 
  
    TES.area1_WAloverlap_mid=48e-6*10e-6;%[m^2]
    TES.area1_WAloverlap_midchk=480e-12;%[m^2]

    TES.area1_WAloverlap_end=258e-12+78.54e-12*2+150e-12+50e-12;%[m^2]
    TES.area1_WAloverlap_endchk=613.719e-12;%[m^2]
else    
    % wafer #9 (PD3) has 20um W/Al overlap
    % for 20um mask
    TES.area1_WAloverlap_mid=48e-6*20e-6;%[m^2]
    TES.area1_WAloverlap_end=1454.87e-12;%[m^2]
end 
TES.vol1_WAloverlap_LEDIT= (2*TES.area1_WAloverlap_end+ 4*TES.area1_WAloverlap_mid)*TES.t; %[m^2]

% Volume of the W only portion of the fin connector
% Since the temperature in the fin connector is lower than the temperature
% in the TES, the effective volume is smaller than the true volume
% This is the efficiency factor for the volume of the fin connector
% contributing to Gep ... we're assuming that this is also the efficiency
% factor for the volume contributing to heat capacity as well.
TES.veff_WFinCon = 0.88;

% Volume of the W/Al overlap portion of the fin connector

% The W/Al portion is completely proximitized ... it should have a very low
% effective volume
%TES.veff_WAloverlap = 0.13; %Measured by Caleb
TES.veff_WAloverlap = 0.45; %We tweak this to match measured bias power (4pW)

if lgc_build
    TES.vol1 =  TES.vol1_TES + TES.veff_WFinCon.*TES.vol1_WFinCon_LEDIT+TES.veff_WAloverlap.*TES.vol1_WAloverlap_LEDIT; %[m^3]
else
    TES.vol1 =  TES.vol1_TES + TES.veff_WFinCon.*TES.vol1_WFinCon+TES.veff_WAloverlap.*TES.vol1_WAloverlap; %[m^3]
end

%---------------- Calculate Resistance ------------------------------------
% Due to electronics constraints, the normal resistance should be ~300mOhm
Rn_wanted = 300e-3; %[Ohm]
%Rn_wanted = 150e-3; %[Ohm]

% while the resistance of a single TES is
Rn_1tes = TES.rho .* TES.l / (TES.w*TES.t); %[Ohm]

% Therefore, the number of TES which can be put in parallel is
%TES.n = ceil(Rn_1tes/Rn_wanted);
%TES.n =1030; % -> let's force this number
%TES.n= 1185;
TES.n= 487+488+28+28;% SZ & MCP:  pulled this from the mask directly!

% Now let's recalculate the Rn value, just a bit more precisely
TES.Rn = Rn_1tes/TES.n; %[Ohm]

%--------------- W Volume -------------------------------------------------
%TES volume
TES.vol= TES.n * TES.vol1; %[m^3]

%----- TES: Operating Point -----
%TES.fOp = 1/3;
%TES.fOp = 1/2;
TES.fOp=0.45;
TES.Ro  = TES.Rn .* TES.fOp; %[Ohm]

%----- Channel Inductance ------
% This is a completely half-ass guess ... I should calculate more carefully
TES.L = 0e-9; %[H]
%TES.L = 25e-9; %[H]


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
wempty = 6e-6;%[m]
wempty_tes = 7.5e-6;%[m]

nhole = 3*QET.nfin;
Ahole = 10e-6*10e-6;%[m]

QET.Afin_empty = QET.nfin*QET.lfin*wempty+ 2*TES.l*wempty_tes+nhole*Ahole; %[m^2]

%QET.Afin = pi* QET.lfin .*(QET.lfin + TES.l/2) - QET.Afin_empty; %[m^2]
QET.Afin = pi* QET.lfin.^2 + 2*QET.lfin*TES.l - QET.Afin_empty; %[m^2]

%we looked at the exact design in L-EDIT (holes need to still be removed)
QET.Afin_LEDIT=28660e-12+28668e-12+28683e-12+28660e-12+28668e-12+28683e-12+nhole*Ahole;%[m^2]

%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.SAactive = Det.nP*TES.n*QET.Afin;%[m^2]

%%%%%%%%%%%%% Passive Aluminum Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rail width
w_rail_main=8e-6;%[m]
w_rail_qet =4e-6;%[m]

%let's calculate the average area per cell:
Acell= Absb.SApattern/(Det.nP*TES.n);%[m^2] -> 1/2 the channels on each side
lcell= sqrt(Acell);%[m]
%let's save
Det.lcell=lcell;%[m]

if lcell > (2*QET.lfin+TES.l)
    %Design isn't close packed
    
    %Passive Al per QET
    Apassive_QET = lcell*w_rail_main + (lcell-(2*QET.lfin+TES.l))*w_rail_qet;%[m^2]
    
else
    % Design is closed packed in TES dimension:
    % there is no vertical rail to the QET:
    ycell = (2*QET.lfin+TES.l); %[m]
    xcell = Acell/ycell;%[m]
    
    Apassive_QET = xcell*w_rail_main;%[m^2]
end    
%Det.SApassive = Apassive_QET*Det.nP*TES.n;

Det.SApassive =   Apassive_QET*Det.nP*TES.n ...                    % TES passive area
                + 2*pi*(Absb.r-Absb.w_safety)*(w_rail_main) ...              % outer ring
                + 2*pi*(Absb.r-Absb.w_safety)*(w_rail_main)/(sqrt(2)) ...    % inner ring
                + 3*(Absb.r-Absb.w_safety)*(1-(sqrt(2)/2))*(w_rail_main) ... % inner vertical rail
                + (Absb.r-Absb.w_safety)*(1+(sqrt(2)/2))*(w_rail_main);      % outer vertical rail

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
