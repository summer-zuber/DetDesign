% Tc plotter -> all devices
Device=[];

%-------- 1cm^3 Si Chip --------
Device.Name     = '1cm3 Si Chip (Design 1 @ UCB) ';
Device.Absb     = 'Si';
Device.Func     = @Chip_1cm3;
Device.Fridge   = fUCB;
Device.elec     = eCDMSII(fUCB);
Device.line_style = '-';
Device.line_color = 'k';

%-------- 1cm^3 Si Chip --------
Device(2).Name     = '1cm3 Si Chip (Design 1 w/ SNOLAB ) ';
Device(2).Absb     = 'Si';
Device(2).Func     = @Chip_1cm3;
Device(2).Fridge   = fSNOLAB;
Device(2).elec     = eSNOLAB(fSNOLAB);
Device(2).line_style = '--';
Device(2).line_color = 'b';


nDev=length(Device);

%Temperature Range to study        
nTc=50;         
Tc=logspace(log10(20e-3),log10(90e-3),nTc); %[K]   

for jDev=1:nDev
    
    Device(jDev).Res_Pt=Tc_ResPt(Tc,Device(jDev).Func,false,Device(jDev).Fridge,Device(jDev).elec,Device(jDev).Absb);
end    

%------ plot -------
leg_str={};

figure(1)
clf(1)
set(1,'position',[25 200 825 600])
hold on
%Scaling Laws
for jDev=1:nDev
    leg_str{jDev}= Device(jDev).Name;
    plot(Tc*1e3,Device(jDev).Res_Pt,[Device(jDev).line_style,Device(jDev).line_color]);
end

% %Measurements
% for jDev=1:nDev
%     nMeas= length(Device(jDev).Meas);
%     for jMeas=1:nMeas
%         leg_str= [leg_str,{Device(jDev).Meas(jMeas).Name}];
%         plot(Device(jDev).Meas(jMeas).Tc*1e3,Device(jDev).Meas(jMeas).sigPt0,[Device(jDev).line_color,Device(jDev).Meas(jMeas).Marker],'Markersize',20);
%     end
% end

hold off
xlim([min(Tc),max(Tc)]*1e3)
set(gca,'yscale','log','xscale','log')
grid on

title('Energy Resolution vs TES Tc')
ylabel('Phonon Energy Resolution (sigma)  [eVt]')
xlabel('TES Tc [mK]')

legend(leg_str,'location','northwest')

set(gca,'xtick',[20:10:100])
xlim([20,100])
ylim([0.1,10])

