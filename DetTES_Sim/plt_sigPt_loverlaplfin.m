function plt_sigPt_loverlaplfin(DetPointer);

if nargin<1
    DetPointer = @(xoverlap,xfin)PhotonDet([],[],[],[],[],xfin,[],xoverlap);
    %DetPointer = @(xoverlap,xfin)     HVZIP([],[],[],[],[],xfin,[],xoverlap);
end

loverlap = [1e-6:2e-6:20e-6]';
n_loverlap=length(loverlap);

%lfin = linspace(100e-6,400e-6,20)';
lfin = [50e-6:10e-6:200e-6]';
n_lfin=length(lfin);

ePcollect  = zeros(n_loverlap,n_lfin);
eQPabsb    = zeros(n_loverlap,n_lfin);
tau_pabsb  = zeros(n_loverlap,n_lfin);
tau_eff    = zeros(n_loverlap,n_lfin);
sigpt      = zeros(n_loverlap,n_lfin);
eEabsb     = zeros(n_loverlap,n_lfin);
Afin       = zeros(n_loverlap,n_lfin);
fAlCovered = zeros(n_loverlap,n_lfin);
for j_lfin=1:n_lfin
    for j_loverlap=1:n_loverlap
        Det=DetPointer(loverlap(j_loverlap),lfin(j_lfin));
        [sigpt(j_loverlap,j_lfin),Det] = ResPt_Estimate(Det);

        ePcollect(j_loverlap,j_lfin) = Det.ePcollect; % Active Al/(Active Al +Passive Al)
        eQPabsb(j_loverlap,j_lfin)   = Det.QET.eQPabsb;
        Afin(j_loverlap,j_lfin)      = Det.QET.Afin;
        tau_pabsb(j_loverlap,j_lfin) = Det.tau_Pabsb;
        tau_eff(j_loverlap,j_lfin)   = Det.TES.tau_etf;
        eEabsb(j_loverlap,j_lfin)    = Det.eEabsb; %-> Total Phonon Energy Collection Taking into account all losses
        fAlCovered(j_loverlap,j_lfin)= Det.fSA_QPabsb; % of total surface area covered by Passive+Active Al
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
title('Energy resolution vs W/Al Overlap Length')
legend(leg_str,'location','northeast')
grid on
ylim([1.5,5])
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
title('% Active Al vs W/Al Overlap Width [um] (for various fin length)')
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
title('% QP Collection Efficiency vs W/AL Overlap (for various fin length)')
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
title('Phonon Collection Bandwidth (for various fin length)')
legend(leg_str,'location','northeast')
grid on

leg_str={};
fig(5)
clf(5)
hold on
for j_lfin=1:n_lfin
    plot(loverlap*1e6,fAlCovered(:,j_lfin)*100,'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('W/Al Overlap Width [um]')
ylabel('% Total Al Coverage [%]')
title('% Total Al Coverage vs W/Al Overlap (for various fin length)')
legend(leg_str,'location','northeast')
grid on
