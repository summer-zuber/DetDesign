function []=RloadJohnsonNoise_Tc()
% let's calculate when the Rload johnson noise is greater than the
% intrinsic TFN noise

TES=[];

TES.nep   = 5;
TES.alpha = 100;
TES.Tb    = [40e-3]';%[K]
TES.To    = linspace(10e-3,100e-3,1e3)';%[K]
TES.Ro    = 90e-3%[Ohm]
TES.Rl    = 32e-3;%[Ohm]
TES.Tl    = 1.4;%[K]
TES.legend='R_{o} = 90mOhm, T_{b}=40mK, \alpha=100';
TES.color = '-k';

TES= NoiseRatio(TES);

TES(2).nep   = 5;
TES(2).alpha = 100;
TES(2).Tb    = [10e-3]';%[K]
TES(2).To    = linspace(10e-3,100e-3,1e3)';%[K]
TES(2).Ro    = 90e-3%[Ohm]
TES(2).Rl    = 32e-3;%[Ohm]
TES(2).Tl    = 1.4;%[K]
TES(2).legend='R_{o} = 90mOhm, T_{b}=10mK, \alpha=100';
TES(2).color ='-r';
TES(2)= NoiseRatio(TES(2));

TES(3).nep   = 5;
TES(3).alpha = 100;
TES(3).Tb    = [10e-3]';%[K]
TES(3).To    = linspace(10e-3,100e-3,1e3)';%[K]
TES(3).Ro    = 180e-3%[Ohm]
TES(3).Rl    = 32e-3;%[Ohm]
TES(3).Tl    = 1.4;%[K]
TES(3).legend='R_{o} = 180mOhm, T_{b}=10mK, \alpha=100';
TES(3).color ='-g';
TES(3)= NoiseRatio(TES(3));

ntes= 3;

figure(102)
clf(102)
hold on
for jplt=1:ntes
    plot(TES(jplt).To(TES(jplt).lgc_valid)*1e3,TES(jplt).LG(TES(jplt).lgc_valid),TES(jplt).color)
end
hold off
grid on
xlabel('T_{c} [mK]')
ylabel('Loop Gain')
legend({TES.legend},'location','southeast')

figure(101)
clf(101)
hold on
for jplt=1:ntes
    plot(TES(jplt).To(TES(jplt).lgc_valid)*1e3,Ftfn(TES(jplt).Tb,TES(jplt).To(TES(jplt).lgc_valid),TES(jplt).nep,true),TES(jplt).color)
end
hold off
grid on
xlabel('T_{c} [mK]')
ylabel('F_{tfn}')
legend({TES.legend},'location','southeast')

figure(99)
clf(99)
hold on
for jplt=1:ntes
    plot(TES(jplt).To(TES(jplt).lgc_valid)*1e3, TES(jplt).rtSpRl_rtSpTFN(TES(jplt).lgc_valid) ,TES(jplt).color)
end
hold off
grid on
xlabel('T_{c} [mK]')
ylabel('R_{load} Johnson Noise / TFN Noise')
legend({TES.legend},'location','southeast')

end

function TES= NoiseRatio(TES)
    TES.LG= TES.alpha./TES.nep.*(1- (TES.Tb./TES.To).^TES.nep);

    %let's only plot those for which the system is stable
    TES.lgc_valid = TES.LG-1 > 0;

    SpRl_SpTFN = (TES.Tl.*TES.Rl)./(TES.To.*TES.Ro).* (1./TES.LG-1).^2 .* (1- (TES.Tb./TES.To).^TES.nep)./ (Ftfn(TES.Tb,TES.To,TES.nep,true).*TES.nep);

    %let's take the sqrt
    TES.rtSpRl_rtSpTFN= sqrt(SpRl_SpTFN);
end


