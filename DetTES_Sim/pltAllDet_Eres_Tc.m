% Tc plotter -> all devices
Device=[];

%-------- Ge iZIP4 @ Berkeley --------
Device.Name     = 'Ge iZIP4 (UCB)';
Device.Absb     = 'Ge';
Device.Func     = @iZIP4;
Device.Fridge   = fUCB;
Device.elec     = eCDMSII(fUCB);
Device.line_style = '-';
Device.line_color = 'g';

%Measurements:
Det= G23R_R475;
Det.Meas.Name= 'G23R Raw';
Det.Meas.Marker='s';
Det.Meas.Tc = Det.TES.Tc;
Device.Meas= Det.Meas;

Det= G23R_R475corr;
Det.Meas.Name= 'G23R Corrected';
Det.Meas.Marker='o';
Det.Meas.Tc = Det.TES.Tc;
Device.Meas(2)= Det.Meas;

 %-------- Ge iZIP4 @ SNOLAB --------
 Device(2).Name = 'Ge iZIP4 (Ideal)';
 Device(2).Absb = 'Ge';
 Device(2).Func = @iZIP4;
 Device(2).Fridge = fSNOLAB;
 Device(2).elec = eSNOLAB;
 Device(2).line_style = '--g';
  
% %-------- SI iZIP4 @ Berkeley --------
Device(3).Name = 'Si iZIP4 (UCB)';
Device(3).Absb = 'Si';
Device(3).Func = @iZIP4;
Device(3).Fridge = fUCB;
Device(3).elec = eCDMSII(fUCB);
Device(3).line_style = '-';
Device(3).line_color = 'k';

%Measurements:
Det= S12C_R466;
Det.Meas.Name= 'S12C Raw';
Det.Meas.Marker='x';
Det.Meas.Tc = Det.TES.Tc;
Device(3).Meas= Det.Meas;

 
 %-------- Si iZIP4 @ SNOLAB --------
 Device(4).Name = 'Si iZIP4 (Ideal)';
 Device(4).Absb = 'Si';
 Device(4).Func = @iZIP4;
 Device(4).Fridge = fSNOLAB;
 Device(4).elec = eSNOLAB;
 Device(4).line_style = '--k';
 
 %-------- Ge HV 3" @ Berkeley --------
 Device(5).Name = 'Ge HV 3"(UCB)';
 Device(5).Absb = 'Ge';
 Device(5).Func = @HV4mm_4ch;
 Device(5).Fridge = fUCB;
 Device(5).elec = eCDMSII(fUCB);
 Device(5).line_style = '-';
 Device(5).line_color = 'y';
  
 %-------- Ge HV 3" @ SNOLAB --------
 Device(6).Name = 'Ge HV 3" (Ideal)';
 Device(6).Absb = 'Ge';
 Device(6).Func = @HV4mm_4ch;
 Device(6).Fridge = fSNOLAB;
 Device(6).elec = eSNOLAB;
 Device(6).line_style = '--';
 Device(6).line_color = 'y';
 

% %-------- Ge SNOLAB iZIP7 --------
% Device(7).Name = 'Ge SNOLAB iZIP ';
% Device(7).Absb = 'Ge';
% Device(7).Func = @iZIP7;
% Device(7).Fridge = fSNOLAB;
% Device(7).elec = eSNOLAB;
% Device(7).line_style = '-c';
% 
% %-------- Si SNOLAB iZIP7 --------
% Device(8).Name = 'Si SNOLAB iZIP ';
% Device(8).Absb = 'Si';
% Device(8).Func = @iZIP7;
% Device(8).Fridge = fSNOLAB;
% Device(8).elec = eSNOLAB;
% Device(8).line_style = '-m';

%-------- Ge SNOLAB HV --------
Device(7).Name = 'Ge SNOLAB HV ';
Device(7).Absb = 'Ge';
Device(7).Func = @HVZIP;
Device(7).Fridge = fSNOLAB;
Device(7).elec = eSNOLAB;
Device(7).line_style = '-';
Device(7).line_color = 'r';
%Let's scale the measurement from iZIP4 design to Ge:
S12C= S12C_R466;
Det=SimplePtScale(S12C,HV2_UMN);
Det.Meas.Name= 'S12C -> Ge HV';
Det.Meas.Marker='x';
Det.Meas.Tc = S12C.TES.Tc;
Device(7).Meas= Det.Meas;

G23R= G23R_R475;
Det=SimplePtScale(G23R,HV2_UMN);
Det.Meas.Name= 'G23R -> Ge HV';
Det.Meas.Marker='s';
Det.Meas.Tc = G23R.TES.Tc;
Device(7).Meas(2)= Det.Meas;

% Det=SimplePtScale(G115f3,HV2_UMN);
% Det.Meas.Name= 'G115f3 -> Ge HV';
% Det.Meas.Marker='d';
% Det.Meas.Tc = G115f3.TES.Tc;
% Device(7).Meas(3)= Det.Meas;

%-------- Si SNOLAB HV --------
Device(8).Name = 'Si SNOLAB HV ';
Device(8).Absb = 'Si';
Device(8).Func = @HVZIP;
Device(8).Fridge = fSNOLAB;
Device(8).elec = eSNOLAB;
Device(8).line_style = '-';
Device(8).line_color = 'b';
%Let's scale the measurement from iZIP4 design to Ge:
S12C= S12C_R466;
Det=SimplePtScale(S12C,SiHV3_UMN);
Det.Meas.Name= 'S12C -> Si HV';
Det.Meas.Marker='x';
Det.Meas.Tc = S12C.TES.Tc;
Device(8).Meas= Det.Meas;

G23R= G23R_R475;
Det=SimplePtScale(G23R,SiHV3_UMN);
Det.Meas.Name= 'G23R -> Si HV';
Det.Meas.Marker='s';
Det.Meas.Tc = G23R.TES.Tc;
Device(8).Meas(2)= Det.Meas;

% Det=SimplePtScale(G115f3,SiHV3_UMN);
% Det.Meas.Name= 'G115f3 -> Si HV';
% Det.Meas.Marker='d';
% Det.Meas.Tc = G115f3.TES.Tc;
% Device(8).Meas(3)= Det.Meas;

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

%Measurements
for jDev=1:nDev
    nMeas= length(Device(jDev).Meas);
    for jMeas=1:nMeas
        leg_str= [leg_str,{Device(jDev).Meas(jMeas).Name}];
        plot(Device(jDev).Meas(jMeas).Tc*1e3,Device(jDev).Meas(jMeas).sigPt0,[Device(jDev).line_color,Device(jDev).Meas(jMeas).Marker],'Markersize',20);
    end
end

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
ylim([5,2e2])

