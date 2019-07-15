% Electro-thermal Oscillation
nRo = 100;
Ro = logspace(log10(20e-3),1,nRo)';%[Ohm]
%Rl = 35e-3; %[Ohm]
Rl = 45e-3; %[Ohm]
L  = 700e-9;%[H]

%gL=20;
gL= 14; % ETF loop gain
beta = 0.2; %[dRdI*I/R]

%W=MaterialProperties('W');
W.gPep_v = 320000000; %[W/K^5/m^3] -> electron-phonon coupling coefficient
W.gC_v = 108;%[J/K^2/m^3] -> metal heat capacity coefficient

%Tc=85e-3;%[K]
Tc=45e-3;
tau0= (W.gC_v .* Tc)./(5 .* W.gPep_v .* Tc.^4);
%tau0=600e-6;%[s]
tau0=78e-6;%[s]

tau_I  = tau0/(1-gL)*ones(nRo,1);%[s]
tau_el = L./(Rl+Ro*(1+beta));%[s]

wavg = 1./tau_I/2+1./tau_el/2;
dw = 1/2.*sqrt((1./tau_el-1./tau_I).^2- 4*Ro/L*gL*(2+beta)/tau0);

wp= wavg+dw; %[rad/s]
wm= wavg-dw; %[rad/s]

taup = 1./wp; %[s]
taum = 1./wm; %[s]

figure(1)
plot(Ro,real(taup)*1e6,'-k')
hold on
plot(Ro,imag(taup)*1e6,'--k')

plot(Ro,real(taum)*1e6,'-b')
plot(Ro,imag(taum)*1e6,'--b')
hold off
xlabel('Ro [Ohm]')
ylabel('time constant [us]')
legend({'Re(\tau_{+})','Im(\tau_{+})','Re(\tau_{-})','Im(\tau_{-})'},'location','best')
grid on
set(gca,'xscale','log')
ylim([-50, 100])
title({'Simple TES poles as a function of Ro',['Loop Gain =', num2str(gL),' beta =',num2str(beta), ' Tc =', num2str(Tc*1e3),'mK']})


figure(2)
plot(Ro,abs(tau_I)*1e6 ,'-r')
hold on
plot(Ro,tau_el*1e6,'-k')
hold off
xlabel('Ro [Ohm]')
ylabel('time constant [us]')
legend({'|\tau_{I}|','\tau_{el}'},'location','best')
set(gca,'xscale','log')
grid on
title({'Simple TES poles as a function of Ro',['Loop Gain =', num2str(gL),' beta =',num2str(beta), ' Tc =', num2str(Tc*1e3),'mK']})



