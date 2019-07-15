function plt_sigPt_loverlaplfin2(DetPointer,Geometry,specs,MatTES);

if nargin<1
    %DetPointer = @(xoverlap,xfin)PhotonDet([],[],[],[],[],xfin,[],xoverlap);
    %DetPointer = @(xoverlap,xfin)     HVZIP([],[],[],[],[],xfin,[],xoverlap);
     DetPointer = @(xGeom,xSpec,xMatTES) Chip_1cm3([],[],xMatTES,[],xGeom,xSpec);
end

if nargin<2
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
    end
    if false
        % Design #1b:  based on the 125umx500um TES square: 1% coverage
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
    end    
    if false
        % Design #1c:  based on the 25umx100um TES square: 1% coverage
        Geometry =[];
        % TES length
        Geometry.ltes = 25e-6;%[m]
        % TES width
        Geometry.wtes = 2e-6;%[m]
        %fin length
        Geometry.lfin = 75e-6; %[m]    
        %fin thickness
        Geometry.hfin = 600e-9;%[m]
        % TES W/Al overlap width
        Geometry.loverlap = 12e-6; %[m]
        % number of TES /channel
        Geometry.ntes= 50;
        % number of fins on a QET
        Geometry.nFin=4;
    end
    if false
        % Design #2
        % R&D upgrades: 1um TES
        % Improved W/Al overlap
        % As W Tc is lowered fraction of W/Al volume that is 
        
        % Design #2:  based on the 25umx50um TES square: 0.2% coverage
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
        Geometry.ntes= 25;
        % number of fins on a QET
        Geometry.nFin=4;
    end
    if false
        % Design #2b
        % R&D upgrades: 1um TES
        % Improved W/Al overlap
        % As W Tc is lowered fraction of W/Al volume that is 
        
        % Design #2:  based on the 25umx50um TES square: 0.2% coverage
        Geometry =[];
        % TES length
        Geometry.ltes = 25e-6;%[m]
        % TES width
        Geometry.wtes = 1e-6;%[m]
        %fin length
        Geometry.lfin = 120e-6; %[m]    
        %fin thickness
        Geometry.hfin = 900e-9;%[m]
        % TES W/Al overlap width
        Geometry.loverlap = 10e-6; %[m]
        % number of TES /channel
        Geometry.ntes= 12;
        % number of fins on a QET
        Geometry.nFin=4;
    end
end

if nargin<3 || isempty(specs)
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

if nargin<4 || isempty(MatTES)
    MatTES= MaterialProperties('W');
    %MatTES.Tc= 20e-3;%[mK]
end

loverlap = [2e-6:1e-6:24e-6]';
n_loverlap=length(loverlap);

%lfin = linspace(100e-6,400e-6,20)';
lfin = [40e-6:20e-6:200e-6]';
n_lfin=length(lfin);

ePcollect  = zeros(n_loverlap,n_lfin);
eQPabsb    = zeros(n_loverlap,n_lfin);
tau_pabsb  = zeros(n_loverlap,n_lfin);
tau_eff    = zeros(n_loverlap,n_lfin);
sigpt      = zeros(n_loverlap,n_lfin);
eEabsb     = zeros(n_loverlap,n_lfin);
Afin       = zeros(n_loverlap,n_lfin);
fSA_APabsb = zeros(n_loverlap,n_lfin);
for j_lfin=1:n_lfin
    for j_loverlap=1:n_loverlap
        Geometry.lfin=lfin(j_lfin);
        Geometry.loverlap=loverlap(j_loverlap);
        Det=DetPointer(Geometry,specs,MatTES);
        
        [sigpt(j_loverlap,j_lfin),Det] = SimulatedNoise_1TES(Det,false);

        ePcollect(j_loverlap,j_lfin) = Det.ePcollect; % Active Al/(Active Al +Passive Al)
        eQPabsb(j_loverlap,j_lfin)   = Det.QET.eQPabsb;
        Afin(j_loverlap,j_lfin)      = Det.QET.A1fin;
        tau_pabsb(j_loverlap,j_lfin) = Det.tau_Pabsb;
        tau_eff(j_loverlap,j_lfin)   = Det.TES.tau_etf;
        eEabsb(j_loverlap,j_lfin)    = Det.eEabsb; %-> Total Phonon Energy Collection Taking into account all losses
        fSA_APabsb(j_loverlap,j_lfin)= Det.fSA_APabsb; % of total surface area covered by Passive+Active Al
    end
end

%/////////////////// Plots ////////////////////////////////////////////////
cmap = jet(n_lfin);

leg_str={};
fig(1)
clf(1)
hold on
for j_lfin=1:n_lfin
    plot(loverlap*1e6,sigpt(:,j_lfin),'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('W/Al Overlap Width [um]')
ylabel('Estimated Baseline Energy Resolution [eV]')
title({'Energy resolution vs W/Al Overlap Length',...
       ['l_{tes}=',num2str(Det.TES.l*1e6),' um Tc = ',num2str(Det.TES.Tc*1e3),' mK']});
legend(leg_str,'location','northeast')
set(gca,'yscale','log')
grid on
ylim([1e-3,1])
xlim([0,20])


leg_str={};
fig(2)
clf(2)
hold on
for j_lfin=1:n_lfin
    plot(loverlap*1e6,ePcollect(:,j_lfin)*100,'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('W/Al Overlap Width [um]')
ylabel('% Active Al [%]')
title({'% Active Al vs W/Al Overlap Width [um] (for various fin length)',...
       ['l_{tes}=',num2str(Det.TES.l*1e6),' um Tc = ',num2str(Det.TES.Tc*1e3),' mK']});
legend(leg_str,'location','northeast')
grid on

leg_str={};
fig(3)
clf(3)
hold on
for j_lfin=1:n_lfin
    plot(loverlap*1e6,eQPabsb(:,j_lfin)*100,'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('W/Al Overlap Width [um]')
ylabel('% QP Collection Efficiency [%]')
title({'% QP Collection Efficiency vs W/AL Overlap (for various fin length)',...
       ['l_{tes}=',num2str(Det.TES.l*1e6),' um Tc = ',num2str(Det.TES.Tc*1e3),' mK']});
legend(leg_str,'location','northeast')
grid on

leg_str={};
fig(4)
clf(4)
hold on
for j_lfin=1:n_lfin
    plot(loverlap*1e6,1./tau_pabsb(:,j_lfin)/(2*pi),'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
for j_lfin=1:n_lfin
    plot(loverlap*1e6,1./tau_eff(:,j_lfin)/(2*pi),'--','color',cmap(j_lfin,:));
end
hold off
xlabel('W/Al Overlap Width [um]')
ylabel('Phonon Collection Bandwidth [hz]')
title({'Phonon Collection Bandwidth (for various fin length)',...
       ['l_{tes}=',num2str(Det.TES.l*1e6),' um Tc = ',num2str(Det.TES.Tc*1e3),' mK']});
legend(leg_str,'location','northeast')
grid on

leg_str={};
fig(5)
clf(5)
hold on
for j_lfin=1:n_lfin
    plot(loverlap*1e6,fSA_APabsb(:,j_lfin)*100,'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('W/Al Overlap Width [um]')
ylabel('% Total Al Coverage [%]')
title({'% Total Al Coverage vs W/Al Overlap (for various fin length)',...
       ['l_{tes}=',num2str(Det.TES.l*1e6),' um Tc = ',num2str(Det.TES.Tc*1e3),' mK']});
legend(leg_str,'location','northeast')
grid on
