function [sigPt_of,Det] = SimulatedNoise_1TES(Det,lgc_plt)
% This function estimates all of the noise sources as a function of frequency
%
%

if nargin <2
    lgc_plt=true;
end    

%Physical Constants
PC= PhysicalConstants;

% Let's calculate the dynamics

if ~isfield(Det.TES,'beta')
    Det.TES.beta=0;
end

if ~isfield(Det.TES,'Qp')
    Det.TES.Qp=0;
end
        
Det = SimpleEquilibrium_1TES(Det,Det.TES.beta,Det.TES.Qp); 

    
% Let's calculate the response functions
Det = Dynamical_Response_1TES(Det,lgc_plt);

omega= Det.Response.omega;
n_omega= length(omega);

if ~isfield(Det.TES,'fSp_xtra');
    Det.TES.fSp_xtra = 0;
    lgc_xtraNoise = false;
else
    lgc_xtraNoise = true;
end

%Squid Noise-----------------------------------------------------
%Let's just assume that the squid noise is frequency invariant 
%Note: this isn't true ... the squid noise clearly get's larger at
%higher frequencies!
Si.Squid = Det.elec.Si_SQUID.^2.* ones(n_omega,1); %[A^2/Hz]
Spt.Squid= Si.Squid ./abs(Det.Response.dIdPt.^2); %[W/Hz]
     
%Johnson Load noise----------------------------------------------
Sv_l = 4*PC.kb* Det.elec.Tl.*Det.elec.Rl; %[V^2/hz]
Si.RlSC = Sv_l./(Det.elec.Rl).^2; %[A^2/Hz]
display(['Superconducting Current Noise is', num2str(round(sqrt(Si.RlSC)*1e12)),'pA/rtHz'])
Si.Rl = abs(Det.Response.dIdV.^2).*Sv_l;%[A^2/Hz]
Spt.Rl = Si.Rl./abs(Det.Response.dIdPt.^2);
    
 %----------------------------------------------------------------
    
 %Johnson TES noise----------------------------------------------
 %Remember that this noise source has both a thermal and an electronic
 %component
  Sv_t= 4*PC.kb*Det.TES.To .* Det.TES.Ro.*(1+Det.TES.beta).^2;%[V^2/Hz]

  Si.Rt  = abs(Det.Response.dIdV- Det.TES.Io.*Det.Response.dIdPt).^2.*Sv_t; %[A^2/Hz]
  Spt.Rt = Si.Rt./abs(Det.Response.dIdPt.^2);%[W^2/Hz]
     
  %---------------------------------------------------------------------
    
  %Phonon cooling noise across TES-Bath Conductance ---------------
    
  %The TFN on a thermal conductance is simply:
  Spt.Gtb = 4*PC.kb*Det.TES.To.^2*Det.TES.Gep.*Ftfn(Det.fridge.T_MC,Det.TES.To,Det.TES.nPep,false).*ones(n_omega,1); %[W^2/rtHz]
  Si.Gtb = abs(Det.Response.dIdPt.^2) .* Spt.Gtb; %[A^2/rtHz]  
  %------------------------------------------------------------------------
  
  %Unexplained Noise that scales as Pt ------------------------------------
  Spt.xtra_Sp = Spt.Gtb * Det.TES.fSp_xtra.^2;
  Si.xtra_Sp  = Si.Gtb * Det.TES.fSp_xtra .^2;
  %------------------------------------------------------------------------
  
%----- Total Noise Terms (Except for 1/f noise) ---------------------------
Si.tot   =  Si.Rl   +Si.Rt   +Si.Gtb   +Si.Squid  + Si.xtra_Sp;
Spt.tot  =  Spt.Rl  +Spt.Rt  +Spt.Gtb  +Spt.Squid + Spt.xtra_Sp;

%-----------Optimum Filter estimators -------------------------------------
%since we may be using non-linear spacing in omega let's find the delta
%omega for every point

domega=zeros(n_omega,1);
domega(2:n_omega-1)= (omega(3:n_omega)-omega(1:n_omega-2))/2;
domega(1) = (omega(2)-omega(1))/2;
domega(n_omega) = (omega(n_omega)-omega(n_omega-1))/2;

sigPt_of_1chan = sqrt(1./sum(domega/(2*pi)*4.*abs(Det.Response.dPtdE).^2./Spt.tot))/PC.qe; %[eV]
sigPt_of = sqrt(Det.nP)*sigPt_of_1chan; %[eV]

display(['Optimum Filter Baseline Resolution for a 100% athermal TES phonon signal in the entire detector is: ',num2str(round(sigPt_of*100)/100),'eVt'])

Det.Noise.omega=omega;
Det.Noise.Si=Si;
Det.Noise.Spt=Spt;
Det.sigPt_of=sigPt_of;


if lgc_plt    

    if lgc_xtraNoise 
        legstr={'Squid','R_{load}','R_{tes}','G: TES-Bath','Extra','Total'};
    else
        legstr={'Squid','R_{load}','R_{tes}','G: TES-Bath','Total'};
    end
        
    gray   = [0.5 ,0.5 ,0.5];
    orange = [0.87,0.49,0];
    peach  = [1.00,0.69,0.39];
    
    %CURRENT NOISE-------------------------------------------------------------
    figure(11)
    plot(omega/(2*pi),sqrt(Si.Squid)*1e12,'-y');
    hold on
    plot(omega/(2*pi),sqrt(Si.Rl)*1e12,'-r');
    plot(omega/(2*pi),sqrt(Si.Rt)*1e12,'-g');
    plot(omega/(2*pi),sqrt(Si.Gtb)*1e12,'-c');
    if lgc_xtraNoise 
        plot(omega/(2*pi),sqrt(Si.xtra_Sp)*1e12,'-b');
    end
    plot(omega/(2*pi),sqrt(Si.tot)*1e12,'-k');
 
    %let's also plot the signal pulse shapes
    yt = Det.Response.dIdPt/abs(Det.Response.dIdPt(1))* 1.25 *(sqrt(Si.tot(1))*1e12);
    plot(omega/(2*pi),abs(yt),'-.','color',gray);
    
    hold off
    xlabel('frequency (hz)')
    ylabel('S_{I} [pA/rtHz')
    set(gca,'xscale','log','yscale','log')
    title('TES Current Noise')
    hleg=legend(legstr,'location','south','fontsize',14);
    set(hleg,'fontsize',20)
    %xlim([1e-2,5e4]);
    %minSi= min(sqrt(Si.tot))*1e12;
    %maxSi= max(sqrt(Si.tot))*1e12;
    %ylim([max(minSi/10,maxSi/5e2),2*maxSi])
    %saveas(jfig+1,['plts/Si_v_',Det.file_str,'.pdf'],'pdf')
    drawnow
    grid on

    %athermal signal NEP ------------------------------------------------------
    figure(12)
    plot(omega/(2*pi),sqrt(Spt.Squid),'-y');
    hold on
    plot(omega/(2*pi),sqrt(Spt.Rl),'-r');
    plot(omega/(2*pi),sqrt(Spt.Rt),'-g');
    plot(omega/(2*pi),sqrt(Spt.Gtb),'-c');
    if lgc_xtraNoise 
        plot(omega/(2*pi),sqrt(Spt.xtra_Sp),'-b');
    end
    plot(omega/(2*pi),sqrt(Spt.tot),'-k');

    %let's also plot the signal pulse shapes
    yt = ones(n_omega,1)*1.25 *(sqrt(Spt.tot(1)));
    plot(omega/(2*pi),abs(yt),'-.','color',gray);
    hold off
    xlabel('frequency (hz)')
    ylabel('S_{p}[W/rtHz]')
    set(gca,'xscale','log','yscale','log')
    title('TES Power Noise')
    hleg=legend(legstr,'location','Southeast','fontsize',14);
    set(hleg,'fontsize',20)
    %xlim([1e-2,5e4]);
    %minSpt= min(sqrt(Spt.tot))*1e15;
    %maxSpt= max(sqrt(Spt.tot))*1e15;
    %ylim([minSpt/100,100*minSpt])
    %saveas(jfig+2,['plts/Spt_v_',Det.file_str,'.pdf'],'pdf')
    grid on
    drawnow
    grid on
end
end

function F=Ftfn(Tl,Th,n,lgcB)
%this function estimates the noise suppression in a thermal conductance
%because the two sides aren't the same temperature

%let's expand out Tl and Th so that they are the same size
nl = length(Tl);
nh = length(Th);

if nl == nh
    %everything is good
elseif nl==1 && nl < nh
    %let's expand out Tl
    Tl = Tl*ones(size(Th));
elseif nh==1 && nh < nl
    %let's expand out Th
    Th = Th*ones(size(Tl));
else
    display('something wrong with Ftfn')
    return 
end    
%let's check to make sure that Tl and Th haven't been switched:
lgc = Tl > Th;
if any(lgc)
    Tx = Th;    
    Th(lgc) = Tl(lgc);
    Tl(lgc) = Tx(lgc);
end

if lgcB
    F = ((Tl./Th).^(n+1)+1)./2;
else    
    F = n/(2*n+1).* ((Tl./Th).^(2*n+1)-1)./((Tl./Th).^(n)-1);
end
end