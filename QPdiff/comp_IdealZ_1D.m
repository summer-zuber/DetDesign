ltes  = 75e-6; %[m]
ltrap = 180e-6;%[m]
lscat = 150e-9;%[m]

%fin variation
nfin=100;
lfin=linspace(.01,3,nfin)'*ltrap;
x=lfin./ltrap;

Afin = 2.*ltes.*lfin; %[m^2]
Amax = 4*lfin(end).^2; %[m^2]

% Code check --------------------------------------------------------------
%let's look at the 0 impedance case ... the two algorithms should be
%identical

fwal= 100; % this is bigger than 1 so it's approximately perfect collection

effabsb_ideal = Effqp_Ideal_rfin(lfin,ltrap);
effabsb_Z = Effqp_Z_rfin(lfin,ltrap,fwal,lscat);

figure(1)
plot(x, effabsb_ideal,'-k')
hold on
plot(x, effabsb_Z,'-b')
hold off
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{qp}')
title('Impedance boundary condition check')

% Let's also check the plot diagnostics
effabsb_Z = Effqp_Z_rfin(180e-6*sqrt(2),ltrap,fwal,lscat,true);
drawnow
effabsb_ideal = Effqp_Ideal_rfin(180e-6*sqrt(2),ltrap,true);
drawnow

% l_fin vs fwal -----------------------------------------------------------
% interface impedance
nf=10;

fwal  = 1./linspace(1,1e4,nf); 

fwal2D= ones(nfin,1)*fwal;
lfin2D= lfin*ones(1,nf);

effabsb_Z = Effqp_Z_rfin(lfin2D,ltrap,fwal2D,lscat);

%------------ Plots -------------------------------------------------------

cmap=colormap(jet(nf));

legstr={};

figure(11)
clf(11)
hold on
for jf=1:nf
    legstr{jf}=['f_{wal}= 1/',num2str(round(1/fwal(jf)))];
    plot(x,effabsb_Z(:,jf)*100,'-','color',cmap(jf,:))
end    
hold off
title('QP collection efficency for 1D Films with W/Al Surface Impedance')
xlabel('l_{fin}/l_{trap}')
ylabel('\epsilon_{collect} %')
ylim([0,100])
legend(legstr,'location','northeast')
