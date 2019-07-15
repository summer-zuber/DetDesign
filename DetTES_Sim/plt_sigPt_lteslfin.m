function plt_sigPt_lteslfin(DetPointer);

if nargin<1
     DetPointer = @(xtes,xfin)PhotonDet([],[],[],[],xtes,xfin,[],[]);
     %DetPointer = @(xtes,xfin)HV4mm_4ch([],[],[],[],xtes,xfin,[],[]);
     %DetPointer = @(xtes,xfin)     HVZIP([],[],[],[],xtes,xfin,[],[]);
end

ltes = [60e-6:20e-6:200e-6]';
n_ltes=length(ltes);

lfin = [100e-6:20e-6:300e-6]';
n_lfin=length(lfin);

ePcollect  = zeros(n_ltes,n_lfin);
eQPabsb    = zeros(n_ltes,n_lfin);
tau_pabsb  = zeros(n_ltes,n_lfin);
tau_eff    = zeros(n_ltes,n_lfin);
sigpt      = zeros(n_ltes,n_lfin);
eEabsb     = zeros(n_ltes,n_lfin);
Afin       = zeros(n_ltes,n_lfin);
fAlCovered = zeros(n_ltes,n_lfin);
for j_lfin=1:n_lfin
    for j_ltes=1:n_ltes
        Det=DetPointer(ltes(j_ltes),lfin(j_lfin));
        [sigpt(j_ltes,j_lfin),Det] = ResPt_Estimate(Det);

        ePcollect(j_ltes,j_lfin) = Det.ePcollect; % Active Al/(Active Al +Passive Al)
        eQPabsb(j_ltes,j_lfin)   = Det.QET.eQPabsb;
        Afin(j_ltes,j_lfin)      = Det.QET.Afin;
        tau_pabsb(j_ltes,j_lfin) = Det.tau_Pabsb;
        tau_eff(j_ltes,j_lfin)   = Det.TES.tau_etf;
        eEabsb(j_ltes,j_lfin)    = Det.eEabsb; %-> Total Phonon Energy Collection Taking into account all losses
        fAlCovered(j_ltes,j_lfin)= Det.fSA_QPabsb; % of total surface area covered by Passive+Active Al
    end
end

%/////////////////// Plots ////////////////////////////////////////////////
cmap = jet(n_lfin);

leg_str={};
fig(1)
clf(1)
hold on
for j_lfin=1:n_lfin
    plot(ltes*1e6,sigpt(:,j_lfin),'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('TES length [um]')
ylabel('Expected Baseline Energy Resolution [eV]')
title('Photon Detector energy resolution vs TES length (for various fin length)')
legend(leg_str,'location','northeast')
grid on





leg_str={};
fig(2)
clf(2)
hold on
for j_lfin=1:n_lfin
    plot(ltes*1e6,ePcollect(:,j_lfin)*100,'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('TES length [um]')
ylabel('% Active Al [%]')
title('% Active Al vs TES length (for various fin length)')
legend(leg_str,'location','northeast')
grid on

leg_str={};
fig(3)
clf(3)
hold on
for j_lfin=1:n_lfin
    plot(ltes*1e6,eQPabsb(:,j_lfin)*100,'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('TES length [um]')
ylabel('% QP Collection Efficiency [%]')
title('% QP Collection Efficiency vs TES length (for various fin length)')
legend(leg_str,'location','northeast')
grid on

leg_str={};
fig(4)
clf(4)
hold on
for j_lfin=1:n_lfin
    plot(ltes*1e6,1./tau_pabsb(:,j_lfin)/(2*pi),'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
%for j_lfin=1:n_lfin
    plot(ltes*1e6,1./tau_eff(:,1)/(2*pi),'--k');
    leg_str=[leg_str,{'Sensor'}]
%end
hold off
xlabel('TES length [um]')
ylabel('Bandwidth [hz]')
title('Phonon Collection and Sensor Bandwidths')
legend(leg_str,'location','northwest')
grid on
ylim([0,2e4])
xlim([50,200])

leg_str={};
fig(5)
clf(5)
hold on
for j_lfin=1:n_lfin
    plot(ltes*1e6,fAlCovered(:,j_lfin)*100,'-','color',cmap(j_lfin,:));
    leg_str{j_lfin}=['l_{fin}=',num2str(round(lfin(j_lfin)*1e6)),'um'];
end
hold off
xlabel('TES length [um]')
ylabel('% Total Al Coverage [%]')
title('% Total Al Coverage vs TES length (for various fin length)')
legend(leg_str,'location','northeast')
grid on

%---------------- SWAP ltes = color, lfin= xaxis --------------------

nplt_ltes = 10;
cmap = jet(nplt_ltes);

indplt_ltes = round(linspace(1,n_ltes,nplt_ltes));

leg_str={};
fig(11)
clf(11)
hold on
for j_ltes= 1:nplt_ltes
    ind = indplt_ltes(j_ltes);
    plot(lfin*1e6,sigpt(ind,:),'-','color',cmap(j_ltes,:));
    leg_str{j_ltes}=['l_{tes}=',num2str(round(ltes(ind)*1e6)),'um'];
end
hold off
xlabel('Al Fin length [um]')
ylabel('Expected Baseline Energy Resolution [eV]')
title('Photon Detector energy resolution vs Fin length (for various TES length)')
legend(leg_str,'location','northeast')
grid on









