% Tc plotter -> all devices
dir_sv = '/Users/mpyle1/CDMS/CDMS_my_papers/UltimatePhononOnly/figs/';

Device=[];

%If you want to vary Al trapping length for new devices ... you can do this
%here:

ldiffQP_old =135e-6; %[m]
ldiffQP_new =600e-6; %[m]

%-------- Ge iZIP4 --------
Device.Name = 'Ge iZIP4';
Device.Absb = 'Ge';
Device.Func = @iZIP4;
Device.line_style = '--r';
Device.ldiffQP =ldiffQP_old;%[m]

%-------- Si iZIP4 --------
Device(2).Name = 'Si iZIP4';
Device(2).Absb = 'Si';
Device(2).Func = @iZIP4;
Device(2).line_style = '--b';
Device(2).ldiffQP =ldiffQP_old;%[m]

%-------- Ge ULZIP --------
Device(3).Name = 'Ge ulZIP';
Device(3).Absb = 'Ge';
Device(3).Func = @ulZIP;
Device(3).line_style = '-r';
Device(3).ldiffQP =ldiffQP_old;%[m]

%-------- Si ULZIP --------
Device(4).Name = 'Si ulZIP';
Device(4).Absb = 'Si';
Device(4).Func = @ulZIP;
Device(4).line_style = '-b';
Device(4).ldiffQP =ldiffQP_old;%[m]

% %-------- Al2O3 ULZIP --------
% Device(5).Name = 'Al2O3 ulZIP';
% Device(5).Absb = 'Al2O3';
% Device(5).Func = @ulZIP;
% Device(5).line_style = '-g';
% Device(5).ldiffQP =ldiffQP_old;%[m]
% 
% %-------- ZnWO4 ULZIP --------
% Device(6).Name = 'ZnWO4 ulZIP';
% Device(6).Absb = 'ZnWO4';
% Device(6).Func = @ulZIP;
% Device(6).line_style = '-c';
% Device(6).ldiffQP =ldiffQP_old;%[m]

nDev=length(Device);

%----------- Experimental Resolutions ----------------------

%-------- Ge G48 Measured --------
MD= G48();

MeasDev=[];

MeasDev.Name = 'G48: Measured ';
MeasDev.line_style = 'o r';
MeasDev.Tc = [MD.Tc];
MeasDev.Res_Pt = [MD.resPt_raw];

MeasDev(2).Name = 'G48: 1/f subtracted';
MeasDev(2).line_style = 'x r';
MeasDev(2).Tc = [MD.Tc];
MeasDev(2).Res_Pt = [MD.resPt_bs];


%-------- Si S12C -------
MD=S12C();

MeasDev(3).Name = 'S12C: Measured';
MeasDev(3).line_style = 'o b';
MeasDev(3).Tc = [MD.Tc];
MeasDev(3).Res_Pt = [MD.resPt_raw];

MeasDev(4).Name = 'S12C: 1/f subtracted';
MeasDev(4).line_style = 'x b';
MeasDev(4).Tc = [MD.Tc];
MeasDev(4).Res_Pt = [MD.resPt_bs];

nMeasDev=length(MeasDev);

% Al properties
Al=MaterialProperties('Al');

%Temperature Range to study        
nTc=50;         
Tc=logspace(log10(10e-3),log10(110e-3),nTc); %[K]   


for jDev=1:nDev
    %let's set the QP length:
    Al.ldiffQP= Device(jDev).ldiffQP;
    
    Device(jDev).Res_Pt=Tc_ResPt(Tc,Device(jDev).Func,false,[],MaterialProperties(Device(jDev).Absb),[],Al);
end    

%------ plot -------
leg_str={};

figure(11)
clf(11)
set(11,'position',[25 200 825 600])
hold on
for jDev=1:nDev
    plot(Tc*1e3,Device(jDev).Res_Pt,Device(jDev).line_style);
    
    %leg_str{jDev}= [Device(jDev).Name,' l_{qp}=',num2str(round(Device(jDev).ldiffQP*1e6)),'um'];
    leg_str{jDev}= [Device(jDev).Name];
end   

for jMD=1:nMeasDev
    plot(MeasDev(jMD).Tc*1e3,MeasDev(jMD).Res_Pt,MeasDev(jMD).line_style,'markersize',10);
 
    leg_str{nDev+jMD}= [MeasDev(jMD).Name];
end   
hold off
xlim([min(Tc),max(Tc)]*1e3)
set(gca,'yscale','log','xscale','log')
grid on

title('Phonon Energy Resolution','FontSize',20,'FontWeight','Bold')
ylabel('Baseline Energy Resolution (sigma)  [eVt]','FontSize',16,'FontWeight','Bold')
xlabel('TES Transition Temperature [mK]','FontSize',16,'FontWeight','Bold')

legend(leg_str,'location','northwest')

print(11,'-depsc',[dir_sv,'EstSigPt_Tc.eps'])

set(gca,'xtick',[20:10:100])
xlim([20,110])
ylim([1e-1,200])

