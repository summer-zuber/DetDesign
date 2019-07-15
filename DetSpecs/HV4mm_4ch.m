function [Det]=HV4mm_4ch(fridge,elec,MatAbsb,MatTES,MatQET,fOp,Qp,lgc_2DOF,Geometry)
% General Detector Object for the light mass WIMP detector (usually with luke neganov gain)
%
% Inputs:
%       1) fridge: fridge object (default fUCB)
%       2) elec:   electronics object (default fUCB)
%       3) MatAbsb: material string/object for the Absorber (default 'Ge')
%       4) MatTES: material string/object for the TES (default: 'W')
%       5) MatQET: material string/object for the QP collection fins (default: 'Al')
%       6) fOp: Operating Point in the transition
%       7) Qp:  Parasitic Power heating of the TES (W)
%       8) lgc_2DOF: treat the TES as if it has internal DOF
%       9) Geometry: Modifications of ltes,lfin,hfin,loverlap
%
% Created MCP (12/3/13)
% 17/04/18: Suhas Ganjam
%      - edited to precisely match design made
%---------------------------------------------------------
if nargin<1 || isempty(fridge)
    %if the fridge is empty, let's use the SNOLAB fridge specs!
    %fridge = fSNOLAB();
    fridge = fUCB(); 
end

if nargin<2 || isempty(elec)
    elec = eCDMSII(fridge); 
end

% Initialize Absorber Object
if nargin<3 || isempty(MatAbsb)
    % if there is no input value, let's assume that the absorber is Ge
    Absb = MaterialProperties('Ge');
elseif iscellstr(MatAbsb)|| ischar(MatAbsb)
    Absb = MaterialProperties(MatAbsb);
else
    % let's inherit all the physical properties from the absorber material
    % to our absorber object
    Absb = MatAbsb;
end  

%Initialize TES Object
if nargin<4 || isempty(MatTES)
    % if there is no input value, let's assume that the TES is W
    TES = MaterialProperties('W');
elseif iscellstr(MatTES)|| ischar(MatTES)
    TES = MaterialProperties(MatTES);
else
    % let's inherit all the physical properties from the TES material
    % to our TES object
    TES = MatTES;
end

%Initialize QET Object
if nargin<5 || isempty(MatQET)
    % if there is no input value, let's assume that the QET fins are made
    % of Al
    QET = MaterialProperties('Al');
elseif iscellstr(MatQET)|| ischar(MatQET)
    QET = MaterialProperties(MatQET);
else
    % let's inherit all the physical properties from fin material
    % to our absorber object
    QET = MatQET;
end

if nargin<6 || isempty(fOp)
    TES.fOp=0.4;
else
    TES.fOp=fOp;
end

if nargin<7 || isempty(Qp)
    TES.Qp = 0;%[W]
else
    TES.Qp = Qp;%[W]
end

if nargin <8 || isempty(lgc_2DOF)
    lgc_2DOF=false;
end

if nargin<9 || isempty(Geometry)
    ltes=200e-6;%[m]
    lfin=240e-6;%[m]
    hfin=600e-9;%[m]
    loverlap=22e-6;%[m]
else
    ltes=Geometry.ltes;
    lfin=Geometry.lfin;
    hfin=Geometry.hfin;
    loverlap=Geometry.loverlap;
end    

%Initialize Detector object
Det=[];
Det.name='G115: 3" HV'; 

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single detector
Det.nP = 4;

%%%%%%%%%%%%%%%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The absorber object inherits all of the properties of the material object
% + ...

% Absorber Dimensions
Absb.h = 4e-3;%[m]
Absb.r = 38.1e-3;%[m]
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
lgc_built = true;
if ~lgc_built
    TES.w = 2.4e-6; %[m]TES width
else
    TES.w = 2.4e-6; %[m]TES width
end

%---------------- Calculate Resistance ------------------------------------
% the normal resistance of a single TES is
Rn_1tes = TES.rho .* TES.l / (TES.w*TES.t); %[Ohm]
    
if ~lgc_built
    % Design Phase:  Let's figure out how many TES we need in an average channel
    
    % Due to electronics constraints, the normal resistance should be ~300mOhm
    Rn_wanted = 350e-3; %[Ohm]

    % Therefore, the number of TES which can be put in parallel is
    TES.n = ceil(Rn_1tes/Rn_wanted);
else
    % Detector Mask has been built and design is fixed
    TES.n = round(3002/4);
end
% Now let's recalculate the Rn value, just a bit more precisely
TES.Rn = Rn_1tes/TES.n; %[Ohm]

%--------------- W Volume -------------------------------------------------

% Volume of a single W TES 
TES.vol1_TES =TES.t*TES.l*TES.w; %[m^3]

% Volume of the W only portion of the fin connectors
TES.nFin = 8; %number of Fins
TES.vol1_WFinCon = TES.nFin*(2.4e-6*10e-6+4e-6*3.4e-6+2*0.75e-6*0.75e-6)*TES.t;%[m^3] non overlapping fin region + neck region connecting fin to TES

% Since the temperature in the fin connector is lower than the temperature
% in the TES, the effective volume is smaller than the true volume
% This is the efficiency factor for the volume of the fin connector
% contributing to Gep ... we're assuming that this is also the efficiency
% factor for the volume contributing to heat capacity as well.
TES.veff_WFinCon =0.88; % From TES Chip Studies

% Volume of the W/Al overlap portion of the fin connector
TES.vol1_WAloverlap=TES.nFin*(2.0e-6*10e-6+ pi*20e-6^2/2)*TES.t;%[m^3] non overlapping fin region + neck region connecting fin to TES

% The W/Al portion is completely proximitized ... it should have a very low
% effective volume
TES.veff_WAloverlap = 0.13; % From TES Chip Studies
%TES.veff_WAloverlap = 0;

%Total effective W volume of a single QET
TES.vol1 = TES.vol1_TES + TES.veff_WFinCon.*TES.vol1_WFinCon+TES.veff_WAloverlap.*TES.vol1_WAloverlap; %[m^3]

% Total effective W volume of a channel
TES.vol  = TES.n*TES.vol1; %[m^3]

%----- TES: Operating Point -----
TES.Ro  = TES.Rn .* TES.fOp; %[Ohm]

%----- Channel Inductance ------
% This is a completely half-ass guess ... I should calculate more carefully
TES.L = 25e-9; %[H]


% %--------------- TES Internal DOF -----------------------------------------
% if lgc_2DOF
%     % for right now let's just put this to be the same as G1b
%     Gtb = TES.nPep .* (TES.gPep_v .* TES.volTES).* TES.Tc.^(TES.nPep-1);%[W/K]
%     G1b = TES.nPep .* (TES.gPep_v .* TES.veff_FinCon.*TES.volFinCon).* TES.Tc.^(TES.nPep-1);%[W/K]
%     
%     TES.nt1=1;
%     TES.Kt1= Gt1/TES.nt1./TES.Tc.^(TES.nt1-1)*2;
% end

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

if ~lgc_built
    % number of fins in a QET
    if QET.lfin > 100e-6
        QET.nfin = 8; %# of fins in a QET
    else
        QET.nfin = 4;
    end
    % width of space between 2 QET fins
    wempty = 6e-6;%[m]
    % width of space around TES
    wempty_tes = 7.5e-6;%[m]

    nhole = 3*QET.nfin;
    Ahole = 10e-6*10e-6;%[m]

    QET.Afin_empty = QET.nfin*QET.lfin*wempty+ 2*TES.l*wempty_tes+nhole*Ahole; %[m^2]

    %QET.Afin = pi* QET.lfin .*(QET.lfin + TES.l/2) - QET.Afin_empty; %[m^2]
    QET.Afin = pi* QET.lfin.^2 + 2*QET.lfin*TES.l - QET.Afin_empty; %[m^2]
else
    QET.nfin = 8;
	fin_areas = [32703.532e-12,27491.371e-12,26931.618e-12,29661.129e-12,32703.513e-12,27491.366e-12,26931.585e-12,29661.185e-12]; %[m^2] each fin has a different area
    QET.Afin = sum(fin_areas);
end

%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.SAactive = Det.nP*TES.n*QET.Afin;%[m^2]

%%%%%%%%%%%%% Passive Aluminum Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backside charge electrode
fSA_charge = .05;

% rail width
w_rail_channel = 8e-6;%[m]
w_rail_main    = 8e-6;%[m]
w_rail_qet     = 4e-6;%[m]

% Area Pads
A_BondPad= 300e-6*700e-6;%[m]
nBondPadTop= 8;
nBondPadBot= 2;

A_SafetyPad = 150e-6^2; %[m]
nSafetyPad  = 11+186+265;

%let's calculate the average area per cell:
Acell= Absb.SApattern/(Det.nP*TES.n);%[m^2] -> all channels on one side
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

if ~lgc_built
	Det.SApassiveTop =   Apassive_QET*Det.nP*TES.n ...                    % TES passive area
						+ 2*pi*(Absb.r-Absb.w_safety)*(w_rail_channel) ...              % outer ring
						+ 2* 2*pi*(Absb.r-Absb.w_safety)*(w_rail_channel)/(sqrt(2)) ...    % inner ring
						+ 8*(Absb.r-Absb.w_safety)*(w_rail_channel) ... % primary rails which make the mercedes design
						+ nBondPadTop*A_BondPad ...
						+ nSafetyPad*A_SafetyPad;
					
	Det.SApassiveBot = Absb.SAface*fSA_charge ...% Bottom side grid%[m^2]
					   + nBondPadBot*A_BondPad; 
else
	Det.SApassiveTop = sum([5177076.116e-12,7405420.785e-12,5626963.311e-12, ...
							6100583.616e-12,5920574.457e-12,6311929.369e-12, ...
							5552807.437e-12,6049754.181e-12]); %Area of all passive Al on top
	
	
	
	Det.SApassiveBot = Absb.SAface*fSA_charge + 2*A_BondPad; %[m^2]
end

Det.SApassive= Det.SApassiveTop+Det.SApassiveBot;

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
%----- Thermal Conductance Coefficient between absorber & Bath -----
Det.Kab = 1.55e-4; %[W/K^4]
%----- Thermal Conductance Coefficient between absorber & Bath
Det.nKab = 4;  

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
