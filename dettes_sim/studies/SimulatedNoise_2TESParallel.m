function [sigPt_of,Det] = SimulatedNoise_2TESParallel(Det,lgc_plt)
% This function estimates all of the noise sources as a function of frequency
%
if nargin <2
    lgc_plt=true;
end    

%Physical Constants
PC= PhysicalConstants;

% Let's calculate the dynamics
Det = SimpleEquilibrium_2TESParallel(Det,Det.TES.beta,Det.TES.Qp); 

% Let's calculate the response functions
Det = Dynamical_Response_2TESParallel(Det,lgc_plt);

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
Sp1.Squid= Si.Squid ./abs(Det.Response.dIdP1.^2); %[W/Hz]
     
%Johnson Load noise----------------------------------------------
Sv_l = 4*PC.kb* Det.elec.Tl.*Det.elec.Rl; %[V^2/hz]
Si.RlSC = Sv_l./(Det.elec.Rl).^2; %[A^2/Hz]
display(['Superconducting Current Noise is', num2str(round(sqrt(Si.RlSC)*1e12)),'pA/rtHz'])
Si.Rl = abs(Det.Response.dIdV.^2).*Sv_l;%[A^2/Hz]
Spt.Rl = Si.Rl./abs(Det.Response.dIdPt.^2);
Sp1.Rl = Si.Rl./abs(Det.Response.dIdP1.^2);
    
 %----------------------------------------------------------------
    
 %Johnson TES noise----------------------------------------------
 %Remember that this noise source has both a thermal and an electronic
 %component
  Sv_t= 4*PC.kb*Det.TES.Tto .* Det.TES.Ro.*(1+Det.TES.beta).^2;%[V^2/Hz]

  Si.Rt  = abs(Det.Response.dIdV- Det.TES.Io.*Det.Response.dIdPt).^2.*Sv_t; %[A^2/Hz]
  Spt.Rt = Si.Rt./abs(Det.Response.dIdPt.^2);%[W^2/Hz]
  Sp1.Rt = Si.Rt./abs(Det.Response.dIdP1.^2);%[W^2/Hz]
        
  %---------------------------------------------------------------------
    
  %Phonon cooling noise across TES-Bath Conductance ---------------
    
  %The TFN on a thermal conductance is simply:
  Spt.Gtb = 4*PC.kb*Det.TES.Tto.^2*Det.TES.Gtb.*Ftfn(Det.fridge.T_MC,Det.TES.Tto,Det.TES.nPep,false).*ones(n_omega,1); %[W^2/Hz]
  Si.Gtb = abs(Det.Response.dIdPt.^2) .* Spt.Gtb; %[A^2/rtHz]  
  Sp1.Gtb = Si.Gtb./abs(Det.Response.dIdP1.^2);%[W^2/Hz];  
  %------------------------------------------------------------------------
  
  %Unexplained Noise that scales as Pt ------------------------------------
  Spt.xtra_Sp = Spt.Gtb * Det.TES.fSp_xtra.^2;
  Sp1.xtra_Sp = Sp1.Gtb * Det.TES.fSp_xtra.^2;
  Si.xtra_Sp  = Si.Gtb * Det.TES.fSp_xtra .^2;
  %------------------------------------------------------------------------
   
  %Phonon cooling noise across 1-Bath Conductance ---------------
  %The TFN on a thermal conductance is simply:
  Sp1.G1b = 4*PC.kb*Det.TES.T1o.^2*Det.TES.G1b.*Ftfn(Det.fridge.T_MC,Det.TES.T1o,Det.TES.nPep,false).*ones(n_omega,1); %[W^2/Hz]
  Si.G1b  = abs(Det.Response.dIdP1.^2) .* Sp1.G1b; %[A^2/rtHz]  
  Spt.G1b = Si.G1b./abs(Det.Response.dIdPt.^2);%[W^2/Hz];  
  %------------------------------------------------------------------------
  
  %Phonon cooling noise across TES-1 Conductance ---------------
  % Here we need to be careful ... energy is conserved, so whatever energy
  % comes out of TES must go inot into C1
  %The TFN on a thermal conductance is simply:
  
  Sq = 4*PC.kb*Det.TES.Tto.^2*Det.TES.G1b.*Ftfn(Det.TES.T1o,Det.TES.Tto,Det.TES.nt1,false).*ones(n_omega,1); %[W^2/Hz]
  
   Si.Gt1 = abs(Det.Response.dIdPt-Det.Response.dIdP1).^2 .* Sq; %[A^2/Hz]  
  Spt.Gt1 = Si.Gt1./abs(Det.Response.dIdPt.^2);%[W^2/Hz];
  Sp1.Gt1 = Si.Gt1./abs(Det.Response.dIdP1.^2);%[W^2/Hz];
  
  %------------------------------------------------------------------------
  
 
%----- Total Noise Terms (Except for 1/f noise) ---------------------------
Si.tot   =  Si.Rl   +Si.Rt   +Si.Gtb   +Si.G1b  + Si.Gt1  +Si.Squid  + Si.xtra_Sp;
Spt.tot  =  Spt.Rl  +Spt.Rt  +Spt.Gtb  +Spt.G1b + Spt.Gt1 +Spt.Squid + Spt.xtra_Sp;
Sp1.tot  =  Sp1.Rl  +Sp1.Rt  +Sp1.Gtb  +Sp1.G1b + Sp1.Gt1 +Sp1.Squid + Sp1.xtra_Sp;

%-----------Optimum Filter estimators -------------------------------------
%since we may be using non-linear spacing in omega let's find the delta
%omega for every point
lgc_pos = omega > 0;

domega=zeros(n_omega,1);
domega(2:n_omega-1)= (omega(3:n_omega)-omega(1:n_omega-2))/2;
domega(1) = (omega(2)-omega(1))/2;
domega(n_omega) = (omega(n_omega)-omega(n_omega-1))/2;

sigPt_of_1chan = sqrt(1./sum(domega(lgc_pos)/(2*pi)*4.*abs(Det.Response.dPtdE(lgc_pos)).^2./Spt.tot(lgc_pos)))/PC.qe; %[eV]
sigPt_of = sqrt(Det.nP)*sigPt_of_1chan; %[eV]

display(['Optimum Filter Baseline Resolution for a 100% athermal TES phonon signal in the entire detector is: ',num2str(round(sigPt_of*100)/100),'eVt'])

Det.Noise.omega=omega;
Det.Noise.Si=Si;
Det.Noise.Spt=Spt;
Det.Noise.Sp1=Sp1;
Det.sigPt_of=sigPt_of;


if lgc_plt    

    if lgc_xtraNoise 
        legstr={'Squid','R_{load}','R_{tes}','G: TES-Bath','G:TES-1','G:1-Bath','Extra','Total'};
    else
        legstr={'Squid','R_{load}','R_{tes}','G: TES-Bath','G:TES-1','G:1-Bath','Total'};
    end
        
    gray   = [0.5 ,0.5 ,0.5];
    orange = [0.87,0.49,0];
    peach  = [1.00,0.69,0.39];
    
    %CURRENT NOISE-------------------------------------------------------------
    figure(11)
    plot(omega(lgc_pos)/(2*pi),sqrt(Si.Squid(lgc_pos))*1e12,'-y');
    hold on
    plot(omega(lgc_pos)/(2*pi),sqrt(Si.Rl(lgc_pos))*1e12,'-r');
    plot(omega(lgc_pos)/(2*pi),sqrt(Si.Rt(lgc_pos))*1e12,'-g');
    plot(omega(lgc_pos)/(2*pi),sqrt(Si.Gtb(lgc_pos))*1e12,'-c');
    plot(omega(lgc_pos)/(2*pi),sqrt(Si.Gt1(lgc_pos))*1e12,'-','color',orange);
    plot(omega(lgc_pos)/(2*pi),sqrt(Si.G1b(lgc_pos))*1e12,'-','color',peach);
    if lgc_xtraNoise 
        plot(omega(lgc_pos)/(2*pi),sqrt(Si.xtra_Sp(lgc_pos))*1e12,'-b');
    end
    plot(omega(lgc_pos)/(2*pi),sqrt(Si.tot(lgc_pos))*1e12,'-k');
 
    %let's also plot the signal pulse shapes
    yt = Det.Response.dIdPt(lgc_pos)/abs(Det.Response.dIdPt(1))* 1.25 *(sqrt(Si.tot(1))*1e12);
    plot(omega(lgc_pos)/(2*pi),abs(yt),'-.','color',gray);
    
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

    %TES signal NEP ------------------------------------------------------
    figure(12)
    plot(omega(lgc_pos)/(2*pi),sqrt(Spt.Squid(lgc_pos)),'-y');
    hold on
    plot(omega(lgc_pos)/(2*pi),sqrt(Spt.Rl(lgc_pos)),'-r');
    plot(omega(lgc_pos)/(2*pi),sqrt(Spt.Rt(lgc_pos)),'-g');
    plot(omega(lgc_pos)/(2*pi),sqrt(Spt.Gtb(lgc_pos)),'-c');
    plot(omega(lgc_pos)/(2*pi),sqrt(Spt.Gt1(lgc_pos)),'-','color',orange);
    plot(omega(lgc_pos)/(2*pi),sqrt(Spt.G1b(lgc_pos)),'-','color',peach);
    
    if lgc_xtraNoise 
        plot(omega(lgc_pos)/(2*pi),sqrt(Spt.xtra_Sp(lgc_pos)),'-b');
    end
    plot(omega(lgc_pos)/(2*pi),sqrt(Spt.tot(lgc_pos)),'-k');

    %let's also plot the signal pulse shapes
    yt = ones(n_omega,1)*1.25 *(sqrt(Spt.tot(1)));
    plot(omega(lgc_pos)/(2*pi),abs(yt(lgc_pos)),'-.','color',gray);
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