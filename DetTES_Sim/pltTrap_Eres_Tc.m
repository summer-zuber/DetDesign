% Tc plotter -> all devices
Device=[];

%If you want to vary Al trapping length for new devices ... you can do this
%here:

ldiffQP= [135e-6,300e-6,600e-6]';

%-------- Ge ulZIP --------
Device.Name = 'Ge ulZIP';
Device.Absb = 'Ge';
Device.Func = @ulZIP;
Device.line_style = '-k';
Device.ldiffQP =ldiffQP(1);%[m]

%-------- Si ulZIP --------
Device(2).Name = 'Si ulZIP';
Device(2).Absb = 'Si';
Device(2).Func = @ulZIP;
Device(2).line_style = '--k';
Device(2).ldiffQP =ldiffQP(1);%[m]

%-------- Ge ulZIP --------
Device(3).Name = 'Ge ulZIP';
Device(3).Absb = 'Ge';
Device(3).Func = @ulZIP;
Device(3).line_style = '-c';
Device(3).ldiffQP =ldiffQP(2);%[m]

%-------- Si ulZIP --------
Device(4).Name = 'Si ulZIP';
Device(4).Absb = 'Si';
Device(4).Func = @ulZIP;
Device(4).line_style = '--c';
Device(4).ldiffQP =ldiffQP(2);%[m]

%-------- Ge ulZIP --------
Device(5).Name = 'Ge ulZIP';
Device(5).Absb = 'Ge';
Device(5).Func = @ulZIP;
Device(5).line_style = '-m';
Device(5).ldiffQP =ldiffQP(3);%[m]

%-------- Si ulZIP --------
Device(6).Name = 'Si ulZIP';
Device(6).Absb = 'Si';
Device(6).Func = @ulZIP;
Device(6).line_style = '--m';
Device(6).ldiffQP =ldiffQP(3);%[m]

nDev=length(Device);

% Al properties
Al=MaterialProperties('Al');

%Temperature Range to study        
nTc=50;         
Tc=logspace(log10(10e-3),log10(100e-3),nTc); %[K]   

leg_str={};

for jDev=1:nDev
    %let's set the QP length:
    Al.ldiffQP= Device(jDev).ldiffQP;
    
    leg_str{jDev}= [Device(jDev).Name,' l_{qp}=',num2str(round(Device(jDev).ldiffQP*1e6)),'um'];
    Device(jDev).Res_Pt=Tc_ResPt(Tc,Device(jDev).Func,false,[],MaterialProperties(Device(jDev).Absb),[],Al);
end    

%------ plot -------
figure(1)
clf(1)
set(1,'position',[25 200 825 600])
hold on
for jDev=1:nDev
    plot(Tc*1e3,Device(jDev).Res_Pt,Device(jDev).line_style);
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
ylim([1e-1,50])

