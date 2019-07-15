function [Det] = Chip_1cm3(fridge,MatAbsb,MatTES,MatQET,Geometry,specs)
% Small Athermal Phonon Chip
%
% Inputs:
%       1) fridge: fridge object (default fSNOLAB)
%       2) MatAbsb: material string/object for the Absorber (default 'Ge')
%       3) MatTES: material string/object for the TES (default: 'W')
%       4) MatQET: material string/object for the QP collection fins (default: 'Al')
%       5) Geometry:
%            ltes
%            wtes
%            lfin
%            hfin
%            loverlap
%            ntes
%       6) Specifications
%           -Quasiparticle trapping efficiency at W/Al interface :
%           eff_absb (deafault = 1.22e-4)
%
% Created MCP 4/14/16)
%---------------------------------------------------------
lgc_diag=true;

if nargin<1 | isempty(fridge)
    %if the fridge is empty, let's use the SNOLAB fridge specs!
    %fridge = fSNOLAB();    
    %fridge = fUCB();
    fridge = fOPD();
end

if strcmpi(fridge.name,'UCB')
    elec= eCDMSII(fridge);
elseif strcmpi(fridge.name,'SNOLAB')
    
    %This device uses low impedance SNOLAB electronics
    elec= eSNOLAB(fridge);
else
    elec= eOPD(fridge);
end

% Initialize Absorber Object
if nargin<2 | isempty(MatAbsb)
    % if there is no input value, let's assume that the absorber is Si
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

if nargin < 5 | isempty(Geometry)
    if true
        % Design #1:  based on the 50umx200um TES square: 1% coverage
        Geometry =[];
        % TES length
        Geometry.ltes = 50e-6;%[m]
        % TES width
        Geometry.wtes = 2e-6;%[m]
        %fin length
        Geometry.lfin = 125e-6; %[m]    
        %fin thickness
        Geometry.hfin = 600e-9;%[m]
        % TES W/Al overlap width
        Geometry.loverlap = 15e-6; %[m]
        % number of TES /channel
        Geometry.ntes= 100;
        % number of fins on a QET
        Geometry.nFin=5;
        % comments:
        %   @ 68mK: 
        %       - 883meV resolution
        %       - tau_Pabsb = 106us
        %       - tau_etf= 18us
        %       - Sensor bandwidth is way larger than absorption bandwidth
        %
        %   @ 40mK: 
        %       - 168meV resolution
        %       - tau_Pabsb = 106us
        %       - tau_etf=96us
        % 
        %
    end
    
    if false
        % Design #1b:  based on the 125umx500um TES square: ?% coverage
        Geometry =[];
        % TES length
        Geometry.ltes = 125e-6;%[m]
        % TES width
        Geometry.wtes = 2e-6;%[m]
        %fin length
        Geometry.lfin = 170e-6; %[m]    
        %fin thickness
        Geometry.hfin = 600e-9;%[m]
        % TES W/Al overlap width
        Geometry.loverlap = 15e-6; %[m]
        % number of TES /channel
        Geometry.ntes= 250;
        % number of fins on a QET
        Geometry.nFin=6;
        % comments:
        %   @ 48mK: 
        %       - 900meV resolution
        %       - tau_Pabsb = 106us
        %       - tau_etf= 18us
        %       - phonon collection can be boosed to 5% without hurting
        %       resolution
        %
        %   @ 45mK: 
        %       - 335meV resolution
        %       - tau_etf=103us
        % 
    end
    
    if false
        % Design #1c:  based on the 25umx100um TES square: ?% coverage
        Geometry =[];
        % TES length
        Geometry.ltes = 25e-6;%[m]
        % TES width
        Geometry.wtes = 2e-6;%[m]
        %fin length
        Geometry.lfin = 100e-6; %[m]    
        %fin thickness
        Geometry.hfin = 600e-9;%[m]
        % TES W/Al overlap width
        Geometry.loverlap = 15e-6; %[m]
        % number of TES /channel
        Geometry.ntes= 50;
        % number of fins on a QET
        Geometry.nFin=4;
        % comments:
        %   @ 48mK: 
        %       - 900meV resolution
        %       - tau_Pabsb = 106us
        %       - tau_etf= 18us
        %       - phonon collection can be boosed to 5% without hurting
        %       resolution
        %
        %   @ 45mK: 
        %       - 335meV resolution
        %       - tau_etf=103us
        % 
    end
    
    if false
        %   Upgraded TES width ->  1um
       
        % Design #2: 25umx50um TES square: 0.8% coverage
        Geometry =[];
        % TES length
        Geometry.ltes = 25e-6;%[m]
        % TES width
        Geometry.wtes = 1e-6;%[m]
        %fin length
        Geometry.lfin = 140e-6; %[m]    
        %fin thickness
        Geometry.hfin = 600e-9;%[m]
        % TES W/Al overlap width
        Geometry.loverlap = 15e-6; %[m]
        % number of TES /channel
        Geometry.ntes= 25;
        % number of fins on a QET
        Geometry.nFin=4;
        % comments:
        %   @ 48mK: 
        %       - 900meV resolution
        %       - tau_Pabsb = 106us
        %       - tau_etf= 18us
        %       - phonon collection can be boosed to 5% without hurting
        %       resolution
        %
        %   @ 45mK: 
        %       - 335meV resolution
        %       - tau_etf=103us
        % 
    end 
    if false
        %   Upgraded TES width ->  1um
       
        % Design #3: 25umx50um TES square: 0.8% coverage
        Geometry =[];
        % TES length
        Geometry.ltes = 25e-6;%[m]
        % TES width
        Geometry.wtes = 1e-6;%[m]
        %fin length
        Geometry.lfin = 140e-6; %[m]    
        %fin thickness
        Geometry.hfin = 900e-9;%[m]
        % TES W/Al overlap width
        Geometry.loverlap = 12e-6; %[m]
        % number of TES /channel
        Geometry.ntes= 12;
        % number of fins on a QET
        Geometry.nFin=4;
    end  
end

if nargin<6 || isempty(specs)
        specs=[];
        %Quasiparticle trapping efficiency
        specs.eff_absb= 1.22e-4;
        %specs.eff_absb = 1e-3;
        
        % Since the temperature in the fin connector is lower than the temperature
        % in the TES, the effective volume is smaller than the true volume
        % This is the efficiency factor for the volume of the fin connector
        % contributing to Gep ... we're assuming that this is also the efficiency
        % factor for the volume contributing to heat capacity as well.

        specs.veff_WAloverlap = 0.13;
        %specs.veff_WAloverlap = 0.05;
        
        % Athermal Phonon Thermalization Probability
        specs.nsurf_pkill   = 400; %[s]
        %specs.nsurf_pkill   = 4000; %[s]
        
end
    
%Initialize Detector object
Det=[];
Det.name='HV Chip';

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single detector
Det.nP = 1;

%%%%%%%%%%%%%%%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The absorber object inherits all of the properties of the material object
% + ...

% Absorber Dimensions
lgc_small=false;
if lgc_small
    Absb.h = 7e-3;%[m]
    Absb.l = 7e-3;%[m]
else
    Absb.h = 10e-3;%[m]
    Absb.l = 10e-3;%[m]
end

% TES patterning isn't possible very near the outer edge
Absb.w_safety = 0e-3; %[m]

Absb.vol  = Absb.l^2 .* Absb.h;  %[m^3]
Absb.SAface = Absb.l^2; %[m^2]

Absb.SApattern = Absb.SAface;%[m^2]
Absb.SA   =  2*Absb.l^2+4*Absb.l.*Absb.h;%[m^2]

% Absorber Mass
Absb.mass = Absb.rhoD_m .* Absb.vol;

% Let's calculate the ballistic scattering length:
Absb.lscat = ScatteringLength_Ballistic(Absb.vol,Absb.SA); %[m]

%%%%%%%%%%%%%%%% TES specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TES.t = 40e-9; %[m] TES thickness
TES.l = Geometry.ltes;%[m]   TES length
TES.w = Geometry.wtes; %[m]TES width

TES.w_Fincon=3.5e-6*0.7;%[m] TES w only finconnector width

% number of TES in a channel 
TES.n= Geometry.ntes;
%---------------- Calculate Resistance ------------------------------------

% single TES resistance
TES.Rn1 = TES.rho .* TES.l / (TES.w*TES.t); %[Ohm]

% Total channel resistance
TES.Rn = TES.Rn1/TES.n; %[Ohm]

%--------------- W Volume -------------------------------------------------
% Volume of a single W TES 
TES.vol1_TES =TES.t*TES.l*TES.w; %[m^3]

% Volume of the W only portion of the fin connectors
TES.nFin = Geometry.nFin; %number of Fins
TES.vol1_WFinCon = TES.nFin * (2*TES.w)*(2*TES.w).*TES.t+ ...  %stub volume
                   2* (TES.l/2) * 2.5e-6 * TES.t; % 1/2 of the end is covered on each side
% if lgc_diag               
%     TES.vol1_WFinCon=0;               
% end

% Since the temperature in the fin connector is lower than the temperature
% in the TES, the effective volume is smaller than the true volume
% This is the efficiency factor for the volume of the fin connector
% contributing to Gep ... we're assuming that this is also the efficiency
% factor for the volume contributing to heat capacity as well.
TES.veff_WFinCon = 0.88; % From TES Chip Studies


% Volume of the W/Al overlap portion of the fin connector
TES.loverlap = Geometry.loverlap; %[m]
%TES.vol1_WAloverlap= 2* (TES.l/2) .* TES.loverlap *TES.t;%[m^3] non overlapping fin region + neck region connecting fin to TES
TES.A1_WAloverlap= TES.nFin*pi/2* TES.loverlap^2;%[m^2] let's use 1/2 circle overlap
TES.vol1_WAloverlap= TES.A1_WAloverlap *TES.t;%[m^3] let's use 1/2 circle overlap

% if lgc_diag
%     TES.vol1_WAloverlap=0;
% end
% The W/Al portion is completely proximitized ... it should have a very low
% effective volume
TES.veff_WAloverlap = specs.veff_WAloverlap; % From TES Chip Studies


%Total effective W volume of a single QET
TES.vol1 = TES.vol1_TES + TES.veff_WFinCon.*TES.vol1_WFinCon+TES.veff_WAloverlap.*TES.vol1_WAloverlap; %[m^3]

% Total effective W volume of a channel
TES.vol  = TES.n*TES.vol1; %[m^3]

%----- TES: Operating Point -----
TES.fOp = 1/3;
%TES.fOp = 1/2;
TES.Ro  = TES.Rn .* TES.fOp; %[Ohm]

%----- Channel Inductance ------
% This is a completely half-ass guess ... I should calculate more carefully
TES.L = 25e-9; %[H]

%%%%%%%%%%%%%%%%%%%%%% QET Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The collection fin object inherits information from its material object
% 2D Diffusion Geometry
QET.lgc_1Ddiff = false;

QET.lfin = Geometry.lfin;%[m]
QET.hfin = Geometry.hfin; %[m]
QET.loverlap = Geometry.loverlap; %[m]

%let's calculate an effective radial overlap assuming that the overlap is a
%1/2 circle with radius loverlap:

% This is the effective radius that we use in our circular symmetery QP
% propagation model
QET.ri = (2*TES.l)/(2*pi);%[m]
QET.roverlap = sqrt((TES.A1_WAloverlap+ pi*QET.ri^2)/pi)-QET.ri;%[m^2] let's use 1/2 circle overlap


QET.eff_absb = specs.eff_absb;
if QET.lgc_1Ddiff
    [QET.eQPabsb,QET.ld,QET.la]= Effqp_1D_moffatt(QET.lfin,QET.hfin,QET.loverlap,QET.eff_absb);
else
    [QET.eQPabsb,QET.ld,QET.la]= Effqp_2D_moffatt(QET.lfin,QET.hfin,QET.roverlap,TES.l,QET.eff_absb);
end

%------ QET Active Area ----- 
% number of fins in a QET
QET.nFin = TES.nFin; %# of fins in a QET

QET.wempty = 4e-6;%[m]
QET.wempty_tes = 8.0e-6;%[m]

fhole = 1e-2; %fraction of Al area which is holes

%QET.Afin = pi* QET.lfin.^2 + 2*QET.lfin*TES.l - QET.Afin_empty; %[m^2]
QET.A1_footprint = pi*QET.lfin*(QET.lfin + TES.l/2);%[m^2]
QET.A1_erase = QET.nFin*QET.wempty*QET.lfin - TES.l*QET.wempty_tes; %[m^2]

QET.A1fin = (QET.A1_footprint-QET.A1_erase)*(1-fhole);%[m^2]

%let's calculate the fraction of the total surface area which is covered by
%QETs
Det.SAactive = Det.nP*TES.n*QET.A1fin;%[m^2]

%%%%%%%%%%%%% Passive Aluminum Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_rail_qet =4e-6;%[m]

%let's close pack all the QETs to limit passive Al

%let's do a ribbon which is 4 QET across:
nQET_width = 4;
n_ribbon= nQET_width + 1;
Det.l_passive= TES.n/nQET_width*(2*QET.lfin); %[m]
Det.SApassive= n_ribbon * w_rail_qet*Det.l_passive;%[m]

Det.SA_Altot= Det.SApassive+Det.SAactive;%[m^2]

% let's caclulate the fraction surface area which is covered with phonon
% absorbing Al
Det.fSA_APabsb = (Det.SApassive+Det.SAactive)./Absb.SA;

% let's calculate the percentage Al surface area which is QET and can thus
% produce signal.
Det.fAl_Active = Det.SAactive/(Det.SApassive+Det.SAactive);

%%%%%%%%%%%%% Ballistic Phonon Collection and Thermalization Time /Efficiency %%%%%%%%%%%%%%%%%%%%

%average athermal phonon collection time in the active and passive Aluminum
Det.tau_pAl = Absb.SA/Det.SA_Altot*Absb.lscat/Absb.vavg_phase;%[s]

% # of surface bouces before thermalization on the bare surface:
Absb.nsurf_pkill   = specs.nsurf_pkill; %[s]
Det.tau_pkill = Absb.nsurf_pkill*Absb.lscat/Absb.vavg_phase;%[s]

Det.tau_Pabsb = 1/(1/Det.tau_pkill + 1/Det.tau_pAl);%[s]
% if lgc_diag
%     Det.tau_Pabsb = Det.tau_Pabsb/10;
% end
Det.w_Pabsb=1/Det.tau_Pabsb; %[rad/s]


% the total fraction of phonon energy which goes into the QP collection fins
% is
Det.ePcollect = Det.tau_pkill/(Det.tau_pkill+Det.tau_pAl) * Det.fAl_Active;

%%%%%%%%%%%%%%% Total Phonon Energy Collection Efficiency %%%%%%%%%%%%%%%%%
 
% The loss mechanisms in our detector are:
% 1) subgap downconversion of athermal phonons in the crystal
% 2) fraction of athermal phonons collected by active metal on the surface of our
% detector instead of collected in passive Al or thermalized (Det.ePcollect) 
% 3) Efficiency of QP production in Al fin (QET.ePQP)
% 4) Efficiency of QP transport to TES (QET.eQPabsb)
% 5) Energy conversion efficiency at W/Al interface
% 6) ?

Det.eE156 = 0.90;

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
