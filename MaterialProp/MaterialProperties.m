function prop = MaterialProperties(MatStr)
% Here we keep track of various material properties!

pc = PhysicalConstants();

prop=[];
if any(strcmpi(MatStr,{'Ge','Germanium'}))
% Germanium Properties
    prop.name='Ge';
    
    prop.lgcElement=true;

    prop.A=72.63; % -> wikipedia
    prop.Z=32;
    
    % Density
    prop.rhoD_m = 5.323e3; %[kg/m^3] -> wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % Ionization Energy
    prop.Eeh = 3.00; % eV; average energy per electron-hole pair (standard definition/measured; assumes electron recoil)
    % bandgap between valence and conduction bands 
    prop.Egap = 0.75;
    % Ionization Fano Factors:
    prop.Fano_ER = 0.13; % from wikipedia: 'fano factor'
    %prop.Fano_ER = 1.0; % from wikipedia: 'fano factor'
    prop.Fano_NR = 1.0;
    
    % The minimum recoil energy required to displace a nucleus in a lattice 
    prop.Er_NRDisp= 15;%[eV] http://scitation.aip.org/content/aip/journal/apl/60/12/10.1063/1.107267
    
    % Debeye Temperature
    prop.Td=374;%[K] from http://link.aps.org/doi/10.1103/PhysRev.124.698
    
    % Phonon Heat Capacity Coefficient
    %prop.gC_mol = 0.38;%[J/kmol/K^4]
    %prop.gC2_v = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
    prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3
    
    %---- Phonon Physics ----
    prop.v_phase =  [3087,3464,5377]';%[m/s]
    
    % Here we average over the density of states assuming spherical symmetry
    % dn/dw = differential density of states with respect to angular
    % frequency = 4*pi* w^2 / v^3
    prop.vavg_phase    = sum(1./prop.v_phase.^2)/sum(1./prop.v_phase.^3);%[m/s] 
    % this is the fraction of longitudinal phonons at a given frequency (isotropic assumption)
    prop.f_long = 1/sum((prop.v_phase(3)./prop.v_phase).^3);
    
    %isotopic scattering rate:  Gamma = a_iso_scat * w^4
    prop.a_iso_scat=  3.67e-41 /(2*pi)^4; %[s^3] from Ge_prop_lattice

    % anharmonic down conversion rate: Gamma = a_ah_dc * w^5
    prop.a_ah_dc = 1.61e-55 /(2*pi)^5; %[s^4] from Ge_prop_lattice
    
    % probability of thermalization when scattering off a surface
    prop.f_surf_dc = 1e-3;
    
    %---- lattice spacing-----
    prop.lconstant = [5.658e-10];%[m]

    %---- Photon Background Ratio (compared to Ge) ---
    prop.fUThK = 1.0; %By Definition!
    
elseif any(strcmpi(MatStr,{'Si','Silicon'}))
%Silicon Properties
    prop.name='Si';
    
    prop.lgcElement=true;

    prop.A=28.085;
    prop.Z=14;
    
    % Density
    prop.rhoD_m = 2.329e3; %[kg/m^3] -> Wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % Ionization Energy
    prop.Eeh = 3.82; % eV; average energy per electron-hole pair (standard definition/measured; assumes electron recoil)
    prop.Egap = 1.4; %eV
    %prop.Fano_ER = .115;  % from wikipedia: 'fano factor'
    prop.Fano_ER = .155;  % from Owens et al, 'On the experimental determination of the Fano factor in Si at soft X-ray wavelengths'
                          %  https://doi.org/10.1016/S0168-9002(02)01178-6
    prop.Fano_NR = 1.0;
    
    % The minimum recoil energy required to displace a nucleus in a lattice 
    prop.Er_NRDisp= 15;%[eV] http://scitation.aip.org/content/aip/journal/apl/60/12/10.1063/1.107267
    
    % Debeye Temperature
    %prop.Td=636;%[K] -> R Keyes paper
    prop.Td=645;%[K] -> wikipedia
    
    %---- lattice spacing-----
    prop.lconstant = [5.431e-10];%[m]
    
    %---- Phonon Physics ----
    % Here are the sound speeds along the [100 direction]
    % http://www.ioffe.ru/SVA/NSM/Semicond/Si/mechanic.html
    prop.v_phase =  [4670,5840,9130]'; %[100 direction]
    
    % Here we average over the density of states assuming spherical symmetry
    % dn/dw = differential density of states with respect to angular
    % frequency = 4*pi* w^2 / v^3
    prop.vavg_phase    = sum(1./prop.v_phase.^2)/sum(1./prop.v_phase.^3);%[m/s] 
    % this is the fraction of longitudinal phonons at a given frequency (isotropic assumption)
    prop.f_long = 1/sum((prop.v_phase(3)./prop.v_phase).^3);
    
     
    %isotopic scattering rate:  Gamma = a_iso_scat * w^4
    % from Tamura, "Quasidiffusive propagation of phonons in silicon: Monte Carlo calculations" , https://doi.org/10.1103/PhysRevB.48.13502
    prop.a_iso_scat=  2.43e6*1e-48 /(2*pi)^4; %[s^3] from Ge_prop_lattice

    % anharmonic down conversion rate: Gamma = a_ah_dc * w^5
    % Tamura, 'Quasidiffusive propagation of phonons in silicon: Monte Carlo calculations', https://doi.org/10.1103/PhysRevB.48.13502
    % We're using the simplified Maris rate calculation
    prop.a_ah_dc = 4.1e4*1e-60 /(2*pi)^5; %[s^4] from Ge_prop_lattice
    
  
    % Phonon Heat Capacity Coefficient
%     prop.gC_mol     =;%[J/kmol/K^4]
%     prop.gC_v     = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
    prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3
    
    %---- Photon Background Ratio (compared to Ge) ---
    prop.fUThK = 2/3; %CDMS-II ratio
    
elseif any(strcmpi(MatStr,{'GaAs'}))
    prop.name='GaAs';
    prop.lgcElement=false;
    
    % Ga has 2 isotopes
    %   69Ga -> 60.11%
    %   71Ga -> 39.89%
    
    % As has 1 isotope
    %   75As -> 100%
    
    % Isotopic scattering will be worse than in Si ... but better than in
    % Ge probably (Ge has 4 different isotopes spanning 70-76)
    prop.Ae = [69.723, 74.921];
    prop.Ze = [31,33];
    prop.Re = [1,1];
    prop.namee={'Ga','As'};
    
    prop.A= sum(prop.Ae.*prop.Re);%AMU
    
    % Density
    prop.rhoD_m = 5.3176e3; %[kg/m^3] -> wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % Debeye Temperature
    prop.Td = 360;%[K] -> wikipedia
    
    % Phonon Heat Capacity Coefficient
%     prop.gC_mol     =;%[J/kmol/K^4]
%     prop.gC_v     = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
      prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3
    
    
    % Ionization Energy
    prop.Eeh = 4.2; % eV; average energy per electron-hole pair (standard definition/measured; assumes electron recoil)
    prop.Egap = 1.4; %eV
    prop.Fano_ER = .12;  % from wikipedia: 'fano factor'
    prop.Fano_NR = 1.0;
    % http://wwwmayr.informatik.tu-muenchen.de/konferenzen/Jass04/courses/4/Tobias%20Eggert/TalkIoffe.pdf
    
    %---- Photon Background Ratio (compared to Ge) ---
    prop.fUThK = 1.0;     
    
    %---- Optical Phonon Frequencies ------------------
    prop.nu_TO= 8.02e12; %[Hz]
    prop.nu_LO= 7.22e12; %[Hz]
    
   
    %---- Phonon Physics ----
    prop.v_phase =  [2480,3350,5240]';%[m/s]
    
    % Here we average over the density of states assuming spherical symmetry
    % dn/dw = differential density of states with respect to angular
    % frequency = 4*pi* w^2 / v^3
    prop.vavg_phase    = sum(1./prop.v_phase.^2)/sum(1./prop.v_phase.^3);%[m/s] 
    % this is the fraction of longitudinal phonons at a given frequency (isotropic assumption)
    prop.f_long = 1/sum((prop.v_phase(3)./prop.v_phase).^3);
    
    %isotopic scattering rate:  Gamma = a_iso_scat * w^4
    prop.a_iso_scat0 = 3.67e-41 /(2*pi)^4; %[s^3] from Ge_prop_lattice
    prop.a_iso_scat1 = 3.67e7/(2*pi*1e12)^4;
    prop.a_iso_scat  = 7.38e6/(2*pi*1e12)^4;%[s^3] "GaAs" Tamura, PRB 31 2574 (1985)

    % anharmonic down conversion rate: Gamma = a_ah_dc * w^5
    prop.a_ah_dc1 = prop.f_long * 1.35e6/(2*pi*1e12)^5; %[s^4] "GaAs-a" Tamura, PRB 31 2574 (1985)
    prop.a_ah_dc = prop.f_long * 7.72e5/(2*pi*1e12)^5; %[s^4] "GaAs-b" Tamura, PRB 31 2574 (1985)
    
    % probability of thermalization when scattering off a surface
    prop.f_surf_dc = 1e-3;
    
elseif any(strcmpi(MatStr,{'TeO2'}))
% Tellerium Oxide -> CUORE
    prop.name='TeO2';
    
    prop.lgcElement=false;
    
    prop.Ae = [127.6, 16];
    prop.Ze = [52,8];
    prop.Re = [1,2];
    prop.namee={'Te','O'};
    
    prop.A= sum(prop.Ae.*prop.Re);%AMU
    
    % Density
    prop.rhoD_m = 6.04e3; %[kg/m^3] -> http://arxiv.org/pdf/1005.3686.pdf
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % Debeye Temperature
    prop.Td=232;%[K] -> http://arxiv.org/pdf/1005.3686.pdf
    
    
    % Phonon Heat Capacity Coefficient
%     prop.gC_mol     =;%[J/kmol/K^4]
%     prop.gC_v     = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
    prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3
    
    %---- Photon Background Ratio (compared to Ge) ---
    prop.fUThK = 1.0; % WHO KNOWS -> THIS IS A PLACE HOLDER
    
    
    
    
elseif any(strcmpi(MatStr,{'CaWO4'}))
% Calcium Tungstate -> CRESST
    prop.name='CaWO4';
    
    prop.lgcElement=false;
 
    prop.Ae = [40.078, 183.84,16];
    prop.Ze = [20,74,8];
    prop.Re = [1,1,4];
    prop.namee={'Ca','W','O'};
    
    prop.A= sum(prop.Ae.*prop.Re);%AMU
    prop.Achk= 287.93; %http://www.guidechem.com/dictionary/en/7790-75-2.html
    
    % Density
    prop.rhoD_m = 6.06e3; %[kg/m^3] -> wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % Debeye Temperature
    prop.Td=250;%[K] -> http://arxiv.org/pdf/0906.3290.pdf
    
    
    % Phonon Heat Capacity Coefficient
%     prop.gC_mol     =;%[J/kmol/K^4]
%     prop.gC_v     = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
    prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3
    
    %---- Photon Background Ratio (compared to Ge) ---
    prop.fUThK = 1.0; % WHO KNOWS -> THIS IS A PLACE HOLDER
    
elseif any(strcmpi(MatStr,{'ZnMO4'}))
% Zinc Molybdate -> Next Gen CUORE
    prop.name='ZnMoO4';
 
    prop.lgcElement=false;
    
    prop.Ae = [65.38, 95.94,16];
    prop.Ze = [30,42,8];
    prop.Re = [1,1,4];
    prop.namee={'Zn','Mo','O'};
    
    prop.A= sum(prop.Ae.*prop.Re);%AMU
    
    % Density
    prop.rhoD_m = 4.317e3; %[kg/m^3] -> http://link.springer.com/article/10.1134%2FS1063774508060266?LI=true#page-1
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % Debeye Temperature
    prop.Td=400;%[K] -> THIS IS JUST A GUESS!!
    
    
    % Phonon Heat Capacity Coefficient
%     prop.gC_mol     =;%[J/kmol/K^4]
%     prop.gC_v     = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
    prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3
    
    %---- Photon Background Ratio (compared to Ge) ---
    prop.fUThK = 1.0; % WHO KNOWS -> THIS IS A PLACE HOLDER
    
elseif any(strcmpi(MatStr,{'ZnWO4'}))
% Zinc Tungsten -> Next Gen CRESST
    prop.name='ZnWO4';
    
    prop.lgcElement=false;
 
    prop.Ae = [65.38, 183.84,16];
    prop.Ze = [30,74,8];
    prop.Re = [1,1,4];
    prop.namee={'Zn','W','O'};
    
    prop.A= sum(prop.Ae.*prop.Re);%AMU
    
    % Density
    prop.rhoD_m = 7.62e3; %[kg/m^3] -> http://www.hilger-crystals.co.uk/prior/mat_znwo4.htm
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % Debeye Temperature
    prop.Td=400;%[K] -> THIS IS JUST A GUESS!!
    
    
    % Phonon Heat Capacity Coefficient
%     prop.gC_mol     =;%[J/kmol/K^4]
%     prop.gC_v     = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
    prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3   
    
    %---- Photon Background Ratio (compared to Ge) ---
    prop.fUThK = 1.0; % WHO KNOWS -> THIS IS A PLACE HOLDER
    
elseif any(strcmpi(MatStr,{'Al2O3','Sapphire'}))
% Sapphire -> CRESST light detector
    prop.name='Al2O3';
    
    prop.lgcElement=false;
    
    %Isotopes:
    % Al:  27Al = 100% 
    % O:  16O = 99.76%
    %     17O =  0.04%
    %     18O =  0.20%
    
    % Very little isotopic scattering
    prop.Ae = [26.981,16];
    prop.Ze = [13,8];
    prop.Re = [2,3];
    prop.namee={'Al','O'};
    
    prop.A= sum(prop.Ae.*prop.Re);%AMU
    
    % Density
    prop.rhoD_m = 3.98e3; %[kg/m^3] -> wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % Debeye Temperature
    prop.Td = 1024;%[K] -> wikipedia
    
    %---- lattice spacing-----
    prop.lconstant = [4.785e-10; 12.991e-10];%  a,c [m]
    
    %---- Phonon Physics ----
    %Let's just use Silicon sound speeds [100 direction]
    % http://www.ioffe.ru/SVA/NSM/Semicond/Si/mechanic.html
    
    % Longitudinal sound speed: https://arxiv.org/pdf/0908.2058.pdf
    % transverse sound speed: 10.1007/BF01325390
    prop.v_phase =  [3/5.5*11000,3/5.5*11000,1*11000]';
    
    % Here we average over the density of states assuming spherical symmetry
    % dn/dw = differential density of states with respect to angular
    % frequency = 4*pi* w^2 / v^3
    prop.vavg_phase    = sum(1./prop.v_phase.^2)/sum(1./prop.v_phase.^3);%[m/s] 
    % this is the fraction of longitudinal phonons at a given frequency (isotropic assumption)
    prop.f_long = 1/sum((prop.v_phase(3)./prop.v_phase).^3);
    
    % Phonon Heat Capacity Coefficient
%     prop.gC_mol     =;%[J/kmol/K^4]
%     prop.gC_v     = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
    prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3
    
    %---- Photon Background Ratio (compared to Ge) ---
    prop.fUThK = 2/3; % Let's use the Si value from CDMS-II
    
elseif any(strcmpi(MatStr,{'W','Tungsten'}))
    %Tungsten
    prop.name='W';
    
    prop.lgcElement=true;
    
    prop.Z = 74;
    prop.A = 183.84;%[AMU] -> wikipedia
    
    % Density
    prop.rhoD_m = 19.25e3; %[kg/m^3] -> wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    % W linear coefficient of heat capacity 
    prop.gC_v = 108; %[J/m^3/K^2] From CRESST LTD7
    
    prop.gC2_mol = 1.01e-3; %[J/gmol/K^2]
    prop.gC2_v = prop.gC2_mol / pc.nA_g * prop.rhoD_n; %[J/m^3/K^2]
    
    %Exponent for linear coefficient of heat capacity;
    prop.nC =1;
    
    % Superconducting heat capacity boost (ratio of heat capacity for superconducting state and normal metal state)
    prop.fCsn =2.43;% let's take into account that in the SC transition the heat capacity is larger

    % electron/phonon thermal coupling coefficient
    prop.gPep_v = 0.32e9; %[W/K^5/m^3] TAREK THESIS p.70
    %prop.gPep_v = 1.8*0.32e9; %[W/K^5/m^3] TAREK THESIS p.70
    
    %From 100x400 40mK TAMU Chip which had 35fW bias power
    % this assumes no parasitic bias power
    prop.gPep_v = 0.22e9; %[W/K^5/m^3]
    
    % Exponent for thermal conductivity between e/h and phonon system
    prop.nPep = 5;
    
    % electronic thermal conductance power law scaling coefficient
    prop.nPee=2; %n-1 is the thermal conductivity scaling coefficient 
    
    % electrical conductivity
    prop.rho_Balzers  = 4.0 * 35e-9; %[Ohm m] -> this is CDMS W from the balzers
    
    % New rho to match the AJA deposition machine at SLAC 9/21/17
    prop.rho_AJA = 1.134*40e-9; %[Ohm m]
    
    % Using G115 to define low Tc TAMU W films
    prop.rho_G115 = 1.5*1.134*40e-9; %[Ohm m]
    
    % MCP 5/7/18: http://titus.stanford.edu/cdms_restricted/detector_physics/HV/ebook/180103/R31_TESChips/eBook/TES_eBook_Run31.html
    %let's use the SLAC TES W Chips (@ 68mK) to measure resistivity:
    prop.rho_SLAC = 0.321*(40e-9*200e-6)/50e-6;%[Ohm m]
    
    %let's use the TAMU TES W Chips (@ 38mK) to measure resistivity:
    prop.rho_TAMU = 0.600*(40e-9*400e-6)/100e-6;%[Ohm m]
    
    prop.rho_Noah  = 1.32e-7; %[Ohm m] -> this matches Noah ... don't know how he got it.
    prop.rho_CRESST = 0.3 * (200e-9*7.5e-3)/5.9e-3; %[Ohm m] -> this is CRESST W
    
    
    prop.rho= prop.rho_TAMU;%[Ohm m]
    %----- Tc ----
    %characteristic Tc will be varied a lot!
    %prop.Tc =68e-3;%[K] -> TES Chip Devices
    %prop.Tc =51e-3; %[K] -> Stanford Si Device
    
    prop.Tc =40e-3;%[K] -> Future Device Guess
    %prop.Tc =15e-3;%[K] -> Future Device Guess
    %prop.Tc=12e-3; %[K] -> Future Device Guess
    
    
    %f_current: this is the normalization factor to BCS theory: i.e. how non perfect the system is in terms of critical current.
    f_current=25;
    prop.dAidAcs = 3.52*sqrt(pc.kb*prop.fCsn*prop.gC_v/pc.hbar/prop.rho)/f_current; %units are totally strange but correct
    
    % Tc width
    %wTc_1090=3.0e-3; %[K] 
    
    %let's scale width with Tc to give us the correct falltime for various
    %devices
    
    % SLAC 68mK W TES test chips
    wTc_1090 =1.4e-3* prop.Tc/68e-3; %[K]
    
    % TAMU 40mK W TES test chips
    % Caleb measured a 66us tau_eff for the 100x400 chip
    wTc_1090 =3.65e-4* prop.Tc/40e-3; %[K]
    
    
    prop.wTc = wTc_1090/2/log(3); %[K]
    
elseif any(strcmpi(MatStr,{'Al','Aluminum'}))
    prop.name='Al';
    
    prop.lgcElement=true;
    
    prop.Z=13;
    prop.A=26.981;
    
    % Efficiency of energy transfer from athermal phonons to quasi-particles 
    prop.ePQP = 0.52;
    
    % QuasiParticle trapping length in our current Al fins
    prop.ldiffQP = 135e-6; %[m]
    %prop.ldiffQP = 300e-6; %[m]
    
    prop.vf = 2.03e6; %[m/s] -> fermi velocity (Ashcroft and Mermin)
    prop.Tc = 1.14; %[K] -> superconducting transition temperature (Ashcroft and Mermin)
    prop.Eg = 1.76*pc.kb.*prop.Tc; %[J] -> cooper pair potential energy
    prop.tau0 =438e-9;%[s]
    
    %tau_al_debeye=3.5e-12; %from kozorezov prb vol 61 num 17 pg 11807
    
    % Debeye Temperature
    prop.Td = 392;%[K] -> Ashcroft and Mermin
    prop.vp_avg = 5100;%[m/s]
    
    % this is the average time to break a cooper pair with a debye
    % frequency phonon:
    prop.tau_ep_debeye= 3.5e-12;%[s] %from kozorezov prb vol 61 num 17 pg 11807
    
    % Resistivity
    prop.RRR = 11; % -> Our iZIPs have an RRR of 11
    %prop.rho_300K = 2.74e-8;%[Ohm m]  http://journals.aps.org/prb/pdf/10.1103/PhysRevB.3.1941
    prop.rho_300K = 27.33e-9;%[Ohm m]  http://journals.jps.jp/doi/pdf/10.1143/JPSJ.66.1253
    prop.rho = prop.rho_300K / prop.RRR;
    
    %number of charge carriers:
    prop.n_e = 18.1e28; %[e/m^3] -> Ashcroft and Mermin
    
    % mass density
    prop.rhoD_m = 2.70e3; %[kg/m^3] -> wikipedia
    
    
    
elseif any(strcmpi(MatStr,{'Ta','Tantalum'}))
    prop.name='Ta';
    
    prop.lgcElement=true;
    
    prop.Z=73;
    prop.A=180.94788;
    
    % Efficiency of energy transfer from athermal phonons to quasi-particles 
    prop.ePQP = 0.52;
    
    prop.Tc = 4.48; %[K] -> superconducting transition temperature (Ashcroft and Mermin)
    prop.Eg = 1.76*pc.kb.*prop.Tc; %[J] -> cooper pair potential energy
    
    % Resistivity
    prop.RRR = 2000; % RRR of 2000 is possible in bulk crystals (http://www.nist.gov/data/PDFfiles/jpcrd258.pdf)
    prop.rho_300K = 131e-9;%[Ohm m]RRR
    prop.rho = prop.rho_300K / prop.RRR;
    
    %number of charge carriers:
    %prop.n_e = 18.1e28; %[e/m^3] -> Ashcroft and Mermin
    
    % mass density
    prop.rhoD_m = 16.69e3; %[kg/m^3] -> wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    %Debeye Temperature
    prop.Td=240;%[K] from wikipedia
    
    % Phonon Heat Capacity Coefficient
    %prop.gC_mol = 0.38;%[J/kmol/K^4]
    %prop.gC_v = prop.gC_mol/prop.A*prop.rhoD_m; % [J/m^3/K^4] = [J/kmol/K^4]*[kmol/kg]*[kg/m^3] 
    prop.gC_v  = 12*pi^4/5*pc.kb*prop.rhoD_n/prop.Td^3; %J/K^4/m^3

    %Speed of Sound
    prop.vP= 3400;%[m/s]
    
elseif any(strcmpi(MatStr,{'Xe','Xenon'}))
%Xenon Properties
    prop.name='Xe';
  
    prop.lgcElement=true;

    prop.A=131.293;
    prop.Z=54;
    
    % Density
    prop.rhoD_m = 3.057e3; %[kg/m^3] -> Wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
elseif any(strcmpi(MatStr,{'He','Helium'}))
%He Properties
    prop.name='He';
  
    prop.lgcElement=true;

    prop.A=4;
    prop.Z=2;
    
    % Density
    prop.rhoD_m = 0.145e3; %[kg/m^3] -> Wikipedia
    prop.rhoD_n = prop.rhoD_m/prop.A*pc.nA_kg; %[#/m^3]
    
    prop.lconstant= (1/prop.rhoD_n/(4/3*pi))^(1/3)*2;%[m]
    
     prop.v_phase =  [239];%[m/s]
else
    display('Pertinent element info not stored')
    display('Please Consider Adding')
    prop = [];
end