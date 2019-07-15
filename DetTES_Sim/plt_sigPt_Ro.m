%let's calculate the noise at a variety of different points along the
%transition

% Detector type
DetObj= @TESrectangles;

% W Properties
W= MaterialProperties('W');
W.Tc= 68.5e-3;%[K]

% Operating Points
% let's go from right above Rl to 0.9Rn;
Det1= DetObj([],[],W);

n_fOp = 5;
fOp=linspace(Det1.elec.Rl./Det1.TES.Rn * 2,0.9,n_fOp)';

%initialize object
Det1= DetObj([],[],W,fOp(1));

Det1.TES.beta=0;
Det1.TES.Qp=0;

[sigPt,Det1]=SimulatedNoise_1TES(Det1,false);
Det=Det1;

%build array
for  jj = 2:n_fOp
    Det1= DetObj([],[],W,fOp(jj));
    Det1.TES.beta=0;
    Det1.TES.Qp=0;
    [sigPt,Det(jj)]=SimulatedNoise_1TES(Det1,false);
end

cmap= jet(n_fOp);
leg_str={};
%let's make some plots!
figure(1)
clf(1)
hold on
for jj=1:n_fOp
    leg_str{jj} = ['R_{0}= ',num2str(round(Det(jj).TES.Ro*1000)),'mOhm'];

    plot(Det(jj).Response.omega/(2*pi),Det(jj).Response.dIdPt ./ Det(jj).Response.dIdV,'-','color',cmap(jj,:))
end    
hold off
xlabel('Frequency [Hz]')
ylabel('dIdP/dIdV [1/A]')
title('dIdP/dIdV for various Ro')
set(gca,'xscale','log','yscale','log')
grid on
legend(leg_str,'location','southwest')

cmap= jet(n_fOp);
leg_str={};
%let's make some plots!
figure(2)
clf(2)
hold on
for jj=1:n_fOp
    leg_str{jj} = ['R_{0}= ',num2str(round(Det(jj).TES.Ro*1000)),'mOhm'];

    plot(Det(jj).Noise.omega/(2*pi),sqrt(Det(jj).Noise.Si.Gtb),'-','color',cmap(jj,:))
end  
for jj=1:n_fOp
    
    plot(Det(jj).Noise.omega/(2*pi),sqrt(Det(jj).Noise.Si.Rl),'--','color',cmap(jj,:))
end  
hold off
xlabel('Frequency [Hz]')
ylabel('Current Noise [A/rtHz]')
title({'Current Noise for various Ro','- TFN , - - R Load'})
set(gca,'xscale','log','yscale','log')
grid on
legend(leg_str,'location','southwest')



