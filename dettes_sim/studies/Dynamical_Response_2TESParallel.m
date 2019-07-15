function Det = Dynamical_Response_2TESParallel(Det,lgc_plt)
% generates the transfer functions for the TES.

% Frequency Range
% let's do the standard fft range so that we can integrate using fft

n_omega=1e6;

nuNyquist= 1e6;%[Hz]
dnu = nuNyquist/n_omega;
nu = [0:n_omega-1]'*dnu;
lgc_nuplt= nu <= nuNyquist/2;
nu(~lgc_nuplt)= nu(~lgc_nuplt)-nuNyquist;

omega = nu*(2*pi);
domega = dnu*(2*pi);

dt= 1/dnu/n_omega;%[s]
t = [0:n_omega-1]'*dt;%[s]
   

Det.Response.omega=omega;
Det.Response.t=t;

if nargin <2
    lgc_plt=true;
end    

%Physical Constants
PC= PhysicalConstants;

%if there is an athermal phonon collection efficiency let's calculate:
lgcAP= isfield(Det,'eEabsb');
if ~lgcAP
    %we want to correctly account for athermal phonon collection
    %efficiency. If this gets used without this ... lets force the numbers
    %to be bad.
    Det.eEabsb=NaN;
    Det.omega_APabsb=1e12; %[1/s] -> this is approximately an infinite bandwidth!
end    
Det.Response.dPtdE= Det.eEabsb./(1+1j.*omega/Det.w_Pabsb);

%if there is no frequency at which 1/f noise
%begins to dominate assume it's not included:
if ~isfield(Det,'lgc_1oF');
    lgc_1oF=false;
else
    lgc_1oF=Det.lgc_1oF;
end

% let's calculate parameters from I.J. Maasilta(1205.5693):

% Parallel Loop Gain
Det.TES.LG = Det.TES.Po .* Det.TES.alpha ./(Det.TES.Gt1_Tt+Det.TES.Gtb)./Det.TES.Tto; %[1]

% Effective Current Biased Positive Feedback Pole
Det.TES.tau_I = Det.TES.Ct./(Det.TES.Gt1_Tt+Det.TES.Gtb)./(1-Det.TES.LG);

% Second dynamics parameter
Det.TES.tau_1 = Det.TES.C1 ./ (Det.TES.Gt1_T1+Det.TES.G1b);

% this is the L/R time constant under the assumption of no pole mixing:L ->
% infty
Det.TES.tau_el  = Det.elec.Lt./(Det.elec.Rl+Det.TES.Ro.*(1+Det.TES.beta));%[s]
Det.TES.w_el = 1./Det.TES.tau_el;%[1/s]

% Calculate dIdV (Using the notation of our fits)
Det.TES.Afit = Det.elec.Rl + Det.TES.Ro.*(1+Det.TES.beta);
Det.TES.Bfit = Det.TES.Ro.*Det.TES.LG./(1-Det.TES.LG).*(2+Det.TES.beta);
Det.TES.Cfit = Det.TES.Gt1_Tt./(Det.TES.Gt1_Tt+Det.TES.Gtb) .* Det.TES.Gt1_T1./(Det.TES.Gt1_T1+Det.TES.G1b)/(1-Det.TES.LG);


% Let's Calculate Jacobian Matrix

%initialize the jacobian matrix
J=zeros(3,3);

%Inductor Equation --------------------------------------------------------
% d(dIdt)dI 
J(1,1)= -1./Det.TES.tau_el;
% d(dIdt)dTt
J(1,2)= - Det.TES.LG.*(Det.TES.Gt1_Tt+Det.TES.Gtb)./(Det.elec.Lt.*Det.TES.Io);
% d(dIdt)dT1
J(1,3)=0;

%TES Energy Dynamical Equation --------------------------------------------
% d(dTtdt)dI 
J(2,1) = Det.TES.Io.*Det.TES.Ro.*(2+Det.TES.beta)./Det.TES.Ct;%[1/s]
% d(dTtdt)dTt 
J(2,2) = -1./Det.TES.tau_I;%[1/s]
% d(dTtdt)dT1 
J(2,3) = Det.TES.Gt1_T1./Det.TES.Ct;%[1/s]

%Heat Capacity 1 Dynamical Equation ---------------------------------------
%d(dT1dt)dI
J(3,1)= 0;
%d(dT1dt)dTt
J(3,2)= Det.TES.Gt1_Tt/Det.TES.C1;
%d(dT1dt)dT1
J(3,3)= -1/Det.TES.tau_1;

% Calculate eigenmodes
%let's find the eigenvalues and match them versus the derived eigenvalues
eig_Jacobian = eig(J); %[1/s]
Det.Response.tau_Jacobian = -1./eig_Jacobian;%[s]

%--------------------------------------------------------------------------

Det.Response.Ztot = Det.TES.Afit.*(1+1j.*omega.*Det.TES.tau_el) +...
                    Det.TES.Bfit./(1+1j.*omega.*Det.TES.tau_I-Det.TES.Cfit./(1+1j.*omega.*Det.TES.tau_1));%[Ohm]
                
Det.Response.Ztes = Det.TES.Ro.*(1+Det.TES.beta) +...
                    Det.TES.Bfit./(1+1j.*omega.*Det.TES.tau_I-Det.TES.Cfit./(1+1j.*omega.*Det.TES.tau_1));%[Ohm]              
 
Det.Response.dIdV = 1./Det.Response.Ztot;


%--------- detM ------------------------
detM =  (  (1+1j.*omega.*Det.TES.tau_el).*(1+1j.*omega.*Det.TES.tau_I).*(1+1j.*omega.*Det.TES.tau_1)-Det.TES.Cfit.*(1+1j.*omega.*Det.TES.tau_el)...
          + Det.TES.Bfit/Det.TES.Afit.*(1+1j.*omega.*Det.TES.tau_1) )...
      ./(Det.TES.tau_el .* Det.TES.tau_I .* Det.TES.tau_1 ); %[1/s^3]

%--------- Calculate dIdVchk -----------
adjM11 = ((1+1j.*omega.*Det.TES.tau_I).*(1+1j.*omega.*Det.TES.tau_1)-Det.TES.Cfit)/(Det.TES.tau_I.*Det.TES.tau_1);
dIdV_chk =  adjM11./detM./Det.elec.Lt;%[Ohm] 

%--------- Calculate dIdPt -------------
adjM12 = - (Det.TES.LG.*(Det.TES.Gt1_Tt+Det.TES.Gtb)./(Det.elec.Lt.*Det.TES.Io)).*((1+1j.*omega.*Det.TES.tau_1)./Det.TES.tau_1);
Det.Response.dIdPt = adjM12./detM./Det.TES.Ct;

dIdPt_chk = -1./(Det.Response.Ztot.*Det.TES.Io).*(Det.Response.Ztes-Det.TES.Ro.*(1+Det.TES.beta))./(Det.TES.Ro.*(2+Det.TES.beta));

%--------- Calculate dIdP1 -------------
adjM13 = -Det.TES.LG./(1-Det.TES.LG).*Det.TES.Gt1_T1 ./(Det.elec.Lt .* Det.TES.Io .* Det.TES.tau_I);
Det.Response.dIdP1 = adjM13./detM./Det.TES.C1;

%--------- Time Domain Calculations ---------------------------------------

% dIdV step function:

dV= (1 - (-1).^round(mod(omega./domega,2)))./(1j.*omega);
dV(1)=0;

Det.Response.dIdV_t= ifft(Det.Response.dIdV .* dV .*dnu*n_omega,'symmetric');

% dIdPt Dirac-Delta pulse shape:
Det.Response.dIdPt_t= ifft(Det.Response.dIdPt.*dnu*n_omega,'symmetric');

% dIdPt Dirac-Delta pulse shape:
Det.Response.dIdP1_t= ifft(Det.Response.dIdP1.*dnu*n_omega,'symmetric');

if lgc_plt
    %------------- dIdV ----------------
    figure(1)
    plot(omega(lgc_nuplt)/(2*pi),abs(Det.Response.dIdV(lgc_nuplt)),'-k')
    hold on
    plot(omega(lgc_nuplt)/(2*pi),abs(dIdV_chk(lgc_nuplt)),'--r')
    %plot(omega/(2*pi), dIdVmid .* ones(n_omega,1),'--c')
    %plot(omega/(2*pi), abs(Det.Response.dIdV0).* ones(n_omega,1),'--c')
    hold off
    xlabel('Frequency [Hz]')
    ylabel('dIdV [1/\Omega]')
    title('Magnitude of dIdV')
    set(gca,'yscale','log','xscale','log')
    grid on

    figure(2)
    plot(omega(lgc_nuplt)/(2*pi),angle(Det.Response.dIdV(lgc_nuplt))*180/pi,'-k')
    hold on
    plot(omega(lgc_nuplt)/(2*pi),angle(dIdV_chk(lgc_nuplt))*180/pi,'--r')
    hold off
    xlabel('Frequency [Hz]')
    ylabel('Phase of dIdV [Deg]')
    title('Phase of dIdV')
    set(gca,'xscale','log')
    grid on
    
     t0 = (t-t(round(n_omega/2)))*1e3;
     lgc_edge = inrange(t0,-2,10);
    
     lgc_end = inrange(t0,4.5,5);
     dIdV_0 = mean(Det.Response.dIdV_t(lgc_end));
     dIdV_t = -(Det.Response.dIdV_t(lgc_edge)-dIdV_0);
     
     figure(3)
     plot(t0(lgc_edge),dIdV_t,'-k')
     %hold on
     %plot(t, Det.Response.dIdV0 .* ones(n_t,1),'--c')
     %plot(t, dIdVmid .* ones(n_t,1),'--c')
     %hold off
     xlabel('time [ms]')
     ylabel('dIdV Step Function [1/Ohm]')
     title('Step function voltage response')
     xlim([-5,5])
     grid on
     
     %let's look at the standard dVdI complex vs real plot
    figure(4)
    set(4,'position',[100 100 800 450])
    plot(real(Det.Response.Ztot(lgc_nuplt)),imag(Det.Response.Ztot(lgc_nuplt)),'-k')
    hold on
    plot(real(Det.Response.Ztes(lgc_nuplt)),imag(Det.Response.Ztes(lgc_nuplt)),'-b')
    hold off
    xlabel('Real(Z) [Ohm]')
    ylabel('Imag(Z) [Ohm]')
    title('Re(Z) vs Im(Z)')
    
    axis equal
    lmax= max(max([abs(real(Det.Response.Ztes)),abs(real(Det.Response.Ztot)),-imag(Det.Response.Ztes)]));
    %let's go to the nearest 10mOhm
    lmax= ceil(lmax*100)/100;
    dx=lmax/5;
    xlim([-lmax,lmax])
    ylim([-lmax,0])
    haxes=gca;
    %let's keep the same grid lines in the x and y direction:
    haxes.XTick= [-lmax:dx:lmax];
    haxes.YTick= [-lmax:dx:0];
    
    legend({'Z_{tot}','Z_{tes}'},'location','north')
    grid on
    
     %------------- dIdPt ----------------
     figure(11)
     plot(omega(lgc_nuplt)/(2*pi),abs(Det.Response.dIdPt(lgc_nuplt)),'-k')
     hold on
     %plot(omega(lgc_nuplt)/(2*pi),abs(dIdPt_chk(lgc_nuplt)),'--r')
     plot(omega(lgc_nuplt)/(2*pi),abs(Det.Response.dIdP1(lgc_nuplt)),'-b')
     hold off
     xlabel('Frequency [Hz]')
     ylabel('|dIdPt| and |dIdP1| [1/V]')
     title('TES Response to Heat Excitations')
     legend({'dIdPt','dIdP1'},'location','northeast')
     set(gca,'yscale','log','xscale','log')
     grid on
     
     figure(12)
     plot(omega(lgc_nuplt)/(2*pi),angle(Det.Response.dIdPt(lgc_nuplt))*180/pi,'-k')
     hold on
     %plot(omega(lgc_nuplt)/(2*pi),angle(dIdPt_chk(lgc_nuplt))*180/pi,'--r')
     plot(omega(lgc_nuplt)/(2*pi),angle(Det.Response.dIdP1(lgc_nuplt))*180/pi,'-b')
     hold off
     xlabel('Frequency [Hz]')
     ylabel('Phase of dIdPt/dIdP1 [Deg]')
     title('TES Response to Heat Stimuli')
     legend({'dIdPt','dIdP1'},'location','northeast')
     set(gca,'xscale','log')
     grid on
     
     figure(13)
     plot(t*1e3,Det.Response.dIdPt_t,'-k')
     hold on
     plot(t*1e3,Det.Response.dIdP1_t,'-b')
     hold off
     xlabel('time [ms]')
     ylabel('dIdPt(t) and dIdP1(t)[1/V]')
     title('TES Response Dirac-Delta Energy Stimuli[1/V]')
     legend({'dIdPt','dIdP1'},'location','northeast')
     xlim([-.1,10])
     grid on
     
     figure(14)
     plot(t*1e3,-Det.Response.dIdPt_t,'-k')
     hold on
     plot(t*1e3,-Det.Response.dIdP1_t,'-b')
     hold off
     xlabel('time [ms]')
     ylabel('dIdPt(t) and dIdP1(t)[1/V]')
     title('TES Response Dirac-Delta Energy Stimuli[1/V]')
     legend({'dIdPt','dIdP1'},'location','northeast')
     xlim([-.1,10])
     set(gca,'yscale','log')
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