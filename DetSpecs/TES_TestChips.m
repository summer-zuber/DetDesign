function [Det]=TES_TestChips(fridge,elec,MatTES,fOp,Qp,lgc_2DOF,Geometry)
% TES Test Chips with Various Fin Connector Designs
%
% Inputs:
%       1) fridge: fridge object (default fUCB)
%       2) elec:   electronics object (default fUCB)
%       3) MatTES: material string/object for the TES (default: 'W')
%       4) fOp: Fractional Operating Point in the transition
%       5) Qp:  Parasitic Power heating of the TES (W)
%       6) lgc_2DOF: treat the TES as if it has internal DOF
%       7) Geometry: 
%           1 (Default): 200um x 800um TES rectangle
%           2:  100um x 400um TES rectangle
%           3:  50um  x 200um TES rectangle
%           4:  25um  x 100um TES rectangle
%           5:  200um TES lines
%           6:  200um TES with 2.5um W only Fin Connector
%           7:  200um TES lines with 2.5um W only Fin Connector and 8um
%               W/Al overlap
%           8:  200um TES lines with 2.5um W only Fin Connector and 20um
%               W/Al overlap
%
% Created MCP (12/3/13)
% 17/04/18: Suhas Ganjam
%      - edited to precisely match design made
%---------------------------------------------------------
if nargin<1 || isempty(fridge)
    %if the fridge is empty, let's use the SNOLAB fridge specs!
    fridge = fSNOLAB();
    %fridge = fUCB(); 
end

if nargin<2 || isempty(elec)
    %elec = eCDMSII(fridge); 
    %elec= eSNOLAB(fridge);
    elec= eSLAC_G115(fridge);
end

%Initialize TES Object
if nargin<3 || isempty(MatTES)
    % if there is no input value, let's assume that the TES is W
    TES = MaterialProperties('W');
elseif iscellstr(MatTES)|| ischar(MatTES)
    TES = MaterialProperties(MatTES);
else
    % let's inherit all the physical properties from the TES material
    % to our TES object
    TES = MatTES;
end

%Fractional Operating Point in the transition
if nargin<4 || isempty(fOp)
    TES.fOp=1/3;
else
    TES.fOp=fOp;
end

% Parasitic Power Heating
if nargin<5 || isempty(Qp)
    TES.Qp = 0;%[W]
else
    TES.Qp = Qp;%[W]
end

if nargin <6 || isempty(lgc_2DOF)
    lgc_2DOF=false;
end

% TES Chip Geometry
if nargin<7 || isempty(Geometry)
%           1 (Default): 200um x 800um TES rectangle
%           2:  100um x 400um TES rectangle
%           3:  50um  x 200um TES rectangle
%           4:  25um  x 100um TES rectangle
%           5:  200um TES lines
%           6:  200um TES with 2.5um W only Fin Connector
%           7:  200um TES lines with 2.5um W only Fin Connector and 8um
%               W/Al overlap
%           8:  200um TES lines with 2.5um W only Fin Connector and 20um
%               W/Al overlap
    Geometry =2;
end    

%Initialize Detector object
Det=[];
 

%%%%%%%%%%%%%%%  Basic DetSpecs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Separately Read Out Phonon Channels on a single detector
Det.nP = 1;

%%%%%%%%%%%%%%%% TES specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Geometry==1
    Det.name='TES Rectangle: 800umx200um';
    
    TES.n = 1; % number of TES in the channel
    TES.nFin = 0; % # of fins on each QET
    
    TES.t = 40e-9; %[m] TES thickness
    TES.l=200e-6;%[m] length of TES (in direction of currect)
    TES.w=800e-6;%[m] width of TES
    
    TES.vol1_WFinCon=0;%[m^3]
    TES.vol1_WAloverlap=0;%[m^3] 
elseif Geometry ==2  
    Det.name='TES Rectangle: 400umx100um';
    
    TES.n = 1; % number of TES in the channel
    TES.nFin = 0; % # of fins on each QET
    
    TES.t = 40e-9; %[m] TES thickness
    TES.l=100e-6;%[m] length of TES (in direction of currect)
    TES.w=400e-6;%[m] width of TES
    
    TES.vol1_WFinCon=0;%[m^3]
    TES.vol1_WAloverlap=0;%[m^3]
elseif Geometry ==3  
    Det.name='TES Rectangle: 200umx50um';
    
    TES.n = 1; % number of TES in the channel
    TES.nFin = 0; % # of fins on each QET
    
    TES.t    = 40e-9; %[m] TES thickness
    TES.l = 500e-6;%[m] length of TES (in direction of currect)
    TES.w = 200e-6;%[m] width of TES
    
    TES.vol1_WFinCon=0;%[m^3]
    TES.vol1_WAloverlap=0;%[m^3]
elseif Geometry ==4  
    Det.name='TES Rectangle: 100umx25um';
    
    TES.n = 1; % number of TES in the channel
    TES.nFin = 0; % # of fins on each QET
    
    TES.t   = 40e-9; %[m] TES thickness
    TES.l   = 25e-6;%[m] length of TES (in direction of currect)
    TES.w   = 100e-6;%[m] width of TES
    
    TES.vol1_WFinCon=0;%[m^3]
    TES.vol1_WAloverlap=0;%[m^3]
elseif Geometry ==5  
    Det.name='200um TES lines';
    
    TES.n = round(800/2.5); % number of TES in the channel
    TES.nFin = 0; % # of fins on each QET
    
    TES.t = 40e-9; %[m] TES thickness
    TES.l = 200e-6;%[m] length of TES (in direction of currect)
    TES.w = 2.5e-6;%[m] width of TES
    
    TES.vol1_WFinCon= 0;%[m^3]
    TES.vol1_WAloverlap=0;%[m^3]
elseif Geometry ==6  
    Det.name='200um TES w/ 2.5um W + 0um W/Al Fin Connector';
    
    TES.n = round(800/2.5); % number of TES in the channel
    TES.nFin = 8; % # of fins on each QET
    
    TES.t = 40e-9; %[m] TES thickness
    TES.l = 200e-6;%[m] length of TES (in direction of currect)
    TES.w = 2.5e-6;%[m] width of TES
    
    TES.vol1_WFinCon= TES.nFin*(  2.5e-6*4e-6 *TES.t ... % [m^3]the neck of the W only fin connector
                               + 2.5e-6*40e-6*TES.t);% [m^3]the long piece W only piece going along the fin
                    
    TES.vol1_WAloverlap=0;%[m^3]
elseif Geometry ==7  
    Det.name='200um TES w/ 2.5um W + 8um W/Al Fin Connector';
    
    TES.n = round(800/2.5); % number of TES in the channel
    TES.nFin = 8; % # of fins on each QET
    
    TES.t = 40e-9; %[m] TES thickness
    TES.l = 200e-6;%[m] length of TES (in direction of currect)
    TES.w = 2.5e-6;%[m] width of TES
    
    TES.vol1_WFinCon = TES.nFin*(  2.5e-6*4e-6 *TES.t ... % [m^3]the neck of the W only fin connector
                                + 2.5e-6*40e-6*TES.t);% [m^3]the long piece W only piece going along the fin
                    
    TES.vol1_WAloverlap = TES.nFin* 40e-6*8e-6*TES.t;%[m^3]
elseif Geometry ==8  
    Det.name='200um TES w/ 2.5um W + 20um W/Al Fin Connector';
    
    TES.n = round(800/2.5); % number of TES in the channel
    TES.nFin = 8; % # of fins on each QET
    
    TES.t = 40e-9; %[m] TES thickness
    TES.l = 200e-6;%[m] length of TES (in direction of currect)
    TES.w = 2.5e-6;%[m] width of TES
    
    TES.vol1_WFinCon = TES.nFin*(  2.5e-6*4e-6 *TES.t ... % [m^3]the neck of the W only fin connector
                                + 2.5e-6*40e-6*TES.t);% [m^3]the long piece W only piece going along the fin
                    
    TES.vol1_WAloverlap = TES.nFin* 40e-6*20e-6*TES.t;%[m^3]
end

%---------------- Calculate Resistance ------------------------------------
% the normal resistance of a single TES is
Rn_1tes = TES.rho .* TES.l / (TES.w*TES.t); %[Ohm]

% Now let's recalculate the Rn value, just a bit more precisely
TES.Rn = Rn_1tes/TES.n; %[Ohm]

%--------------- W Volume -------------------------------------------------
%TES volume

% Volume of a single W TES 
TES.vol1_TES =(TES.t*TES.l*TES.w); %[m^3]

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
TES.veff_WAloverlap = 0.13;

TES.vol =  TES.n*(TES.vol1_TES + TES.veff_WFinCon.*TES.vol1_WFinCon+TES.veff_WAloverlap.*TES.vol1_WAloverlap); %[m^3]

%----- TES: Operating Point -----
TES.Ro  = TES.Rn .* TES.fOp; %[Ohm]

%----- Channel Inductance ------
% This is a completely half-ass guess ... I should calculate more carefully
TES.L = 5e-9; %[H]


%--------------- TES Internal DOF -----------------------------------------
if lgc_2DOF
    % for right now let's just put this to be the same as G1b
    Gtb = TES.nPep .* (TES.gPep_v .* TES.volTES).* TES.Tc.^(TES.nPep-1);%[W/K]
    G1b = TES.nPep .* (TES.gPep_v .* TES.veff_FinCon.*TES.volFinCon).* TES.Tc.^(TES.nPep-1);%[W/K]
    
    TES.nt1=1;
    TES.Kt1= Gt1/TES.nt1./TES.Tc.^(TES.nt1-1)*2;
end

%%%%%%%%%%%%%%%%%%% Thermal Conductance to Bath %%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Thermal Conductance Coefficient between absorber & Bath -----
Det.Gab.K = 1.55e-4; %[W/K^4]
%----- Thermal Conductance Coefficient between absorber & Bath
Det.Gab.n = 4;  

%%%%%%%%%%%%%%%% Electronics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let's sum all the inductances to find the total inductance
elec.Lt = elec.Lsquid+elec.Lp+TES.L; %[H]

%//////////////// Measured Quantities /////////////////////
% Baseline energy resolution
Det.ResPt =[];%[sigma eVt]

%%%%%%%%%%%%%%%% Combine Objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Det.TES=TES;
Det.elec=elec;
Det.fridge=fridge;
