%let's do some checks on these codes

%--- Check 1: vary fin length ---
% let's vary Al fin length, keeping overlap constant, fin height
lfin = [100e-6:10e-6:2e-3]';%[m]
ltes = 200e-6;%[m]
hfin = 600e-9;%[m]
loverlap = 20e-6;%[m]

figure(1)
[eff_1D,ld_1D,la_1D,effabsb_1D] = Effqp_1D_moffatt(lfin,hfin,loverlap);
[eff_2D,ld_2D,la_2D,effabsb_2D] = Effqp_2D_moffatt(lfin,hfin,loverlap,ltes);
plot(lfin*1e6,eff_1D,'-k')
hold on
plot(lfin*1e6,eff_2D,'-r')
hold off
xlabel('fin length [um]')
ylabel('QP collection probability')
title({'QP collection: vary fin length',['h_{fin}=',num2str(round(hfin*1e8)/1e2),'um l_{w/al}=',num2str(round(loverlap*1e8)/1e2),'um l_{tes}=',num2str(round(ltes*1e8)/1e2),'um']})
legend({'1D','2D'},'location','northeast')
grid on
saveas(1,'plt/effQP_lfin.png','png')

%--- Check 2: vary fin height ---
% let's vary Al fin height, keeping overlap constant, fin height
lfin = [200e-6];%[m]
ltes = 200e-6;%[m]
hfin = [100e-9:10e-9:2e-6]';%[m]
loverlap = 20e-6;%[m]

figure(11)
[eff_1D,ld_1D,la_1D,effabsb_1D] = Effqp_1D_moffatt(lfin,hfin,loverlap);
[eff_2D,ld_2D,la_2D,effabsb_2D] = Effqp_2D_moffatt(lfin,hfin,loverlap,ltes);
plot(hfin*1e9,eff_1D,'-k')
hold on
plot(hfin*1e9,eff_2D,'-r')
hold off
xlabel('Al fin thickness [nm]')
ylabel('QP collection probability')
title({'QP collection: vary fin thickness',['l_{fin}=',num2str(round(lfin*1e8)/1e2),'um l_{w/al}=',num2str(round(loverlap*1e8)/1e2),'um l_{tes}=',num2str(round(ltes*1e8)/1e2),'um']})
grid on
legend({'1D','2D'},'location','southeast')
saveas(11,'plt/effQP_hfin.png','png')

figure(12)
plot(hfin*1e9,ld_1D*1e6,'-k')
hold off
xlabel('Al fin thickness [nm]')
ylabel('Diffusion Length')
title({'QP Diffusion Length: vary fin thickness',['l_{fin}=',num2str(round(lfin*1e8)/1e2),'um l_{w/al}=',num2str(round(loverlap*1e8)/1e2),'um l_{tes}=',num2str(round(ltes*1e8)/1e2),'um']})
grid on

%---- Check: vary overlap ----
% let's vary Al fin height, keeping overlap constant, fin height
lfin = 200e-6;%[m]
ltes = 200e-6;%[m]
hfin = 600e-9;%[m]
loverlap = [1e-6:1e-6:100e-6]';%[m]

figure(21)
[eff_1D,ld_1D,la_1D,effabsb_1D] = Effqp_1D_moffatt(lfin,hfin,loverlap);
[eff_2D,ld_2D,la_2D,effabsb_2D] = Effqp_2D_moffatt(lfin,hfin,loverlap,ltes);
plot(loverlap*1e6,eff_1D,'-k')
hold on
plot(loverlap*1e6,eff_2D,'-r')
xlabel('W/Al overlap width[um]')
ylabel('QP collection probability')
title({'QP collection: vary W/Al overlap',['l_{fin}=',num2str(round(lfin*1e8)/1e2),'um h_{fin}=',num2str(round(hfin*1e8)/1e2),'um l_{tes}=',num2str(round(ltes*1e8)/1e2),'um']})
grid on
legend({'1D','2D'},'location','southeast')
saveas(21,'plt/effQP_loverlap.png','png')

%---- Check: vary TES length ----
% let's vary TES and Fin length
lfin = [0e-6:10e-6:1000e-6];%[m]

n_ltes =4;
%ltes = logspace(-4,-1,n_ltes)';%[m]
ltes=[100e-6,300e-6,500e-6,100e-3]';%[m]

lfinM= repmat(lfin,size(ltes));
ltesM= repmat(ltes,size(lfin));

hfin = 600e-9;%[m]
%loverlap = 5e-6;%[m]
loverlap = 20e-6;%[m]

[eff_1D,ld_1D,la_1D,effabsb_1D] = Effqp_1D_moffatt(lfin,hfin,loverlap);
[eff_2D,ld_2D,la_2D,effabsb_2D] = Effqp_2D_moffatt(lfinM,hfin,loverlap,ltesM);

cmap=jet(n_ltes);
leg_str={};

figure(31)
clf(31)
hold on
for j_ltes=1:n_ltes
    leg_str(j_ltes)={['l_{tes}=',num2str(round(ltes(j_ltes)*1e6)),'um']}
    plot(lfinM(j_ltes,:)*1e6,eff_2D(j_ltes,:),'-','color',cmap(j_ltes,:));
end
plot(lfin*1e6,eff_1D,'--k');
leg_str(n_ltes+1)={'1D QP'};
hold off
xlabel('Al fin length [um]')
ylabel('QP collection probability')
title({'2D QP collection: vary l_{tes} & l_{fin}',['h_{fin}=',num2str(round(hfin*1e8)/1e2),'um l_{overlap}=',num2str(round(loverlap*1e8)/1e2),'um']})
grid on
legend(leg_str,'location','northeast')
saveas(31,'plt/effQP_ltes.png','png')



