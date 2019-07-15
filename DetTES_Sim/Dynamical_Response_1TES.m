function Det = Dynamical_Response_1TES(Det,lgc_plt,omega)
% generates the transfer functions for the TES.
%

% Let's calculate the dynamics
%Det=SimpleEquilibrium_1TES(Det,Det.TES.beta,Det.TES.Qp);	

%Frequency Range
if nargin<3
    %n_omega=2.5e5;
    n_omega=1e5;
    %n_omega=1e2;
    omega= logspace(-1,5.5,n_omega)'*(2*pi);
    %omega = linspace(0,1e5,n_omega)'*(2*pi);
else
    n_omega=length(omega);
end
Det.Response.omega=omega;

if nargin <2
    lgc_plt=true;
end    

%plot simple approximations
lgc_pltsimp = true;

%Physical Constants
PC= PhysicalConstants;


%if there is an athermal phonon collection efficiency let's calculate:
lgcAP= isfield(Det,'eEabsb');
if ~lgcAP
    %we want to correctly account for athermal phonon collection
    %efficiency. If this gets used without this ... lets force the numbers
    %to be bad.
    Det.eEabsb=1;
    Det.w_Pabsb=1e12; %[1/s] -> this is approximately an infinite bandwidth!
end    
Det.Response.dPtdE= Det.eEabsb./(1+1j.*omega/Det.w_Pabsb);

%if there is no frequency at which 1/f noise
%begins to dominate assume it's not included:
if ~isfield(Det,'lgc_1oF');
    lgc_1oF=false;
else
    lgc_1oF=Det.lgc_1oF;
end

% Calculate dIdP
Det.Response.dIdPt = - Det.TES.Gep .* Det.TES.LG ./(Det.TES.C*Det.TES.Io*Det.elec.Lt)./...
    ( (1j.*omega + Det.TES.Gep*(1-Det.TES.LG)./Det.TES.C).*(1j.*omega + (Det.elec.Rl+Det.TES.Ro.*(1+Det.TES.beta))/Det.elec.Lt)+Det.TES.LG*Det.TES.Ro/Det.elec.Lt*Det.TES.Gep/Det.TES.C.*(2+Det.TES.beta));

% Calculate dIdV
Det.Response.Ztes = Det.TES.Ro.*(1+Det.TES.beta)+...
     Det.TES.Ro.*Det.TES.LG./(1-Det.TES.LG).*(2+Det.TES.beta)./(1+1j.*omega.*Det.TES.C./Det.TES.Gep./(1-Det.TES.LG));%[Ohm]
 
Det.Response.Ztot = Det.elec.Rl+...
     1j.*omega.*Det.elec.Lt+...
     Det.Response.Ztes;%[Ohm]

 Det.Response.dIdV = 1./Det.Response.Ztot;


% Let's Calculate dIdV step functions
%Using Irwin's notation

% C/G for the simple thermal circuit
Det.TES.tau0 = Det.TES.C ./ Det.TES.Gep; %[s]

% exponetial rise time for the current biased circuit 
Det.TES.tau_I = Det.TES.tau0./(1-Det.TES.LG); %[s]
Det.TES.w_I = 1./Det.TES.tau_I;%[1/s]

% this is the L/R time constant under the assumption of no pole mixing:L ->
% infty
Det.TES.tau_el  = Det.elec.Lt./(Det.elec.Rl+Det.TES.Ro.*(1+Det.TES.beta));%[s]
Det.TES.w_el = 1./Det.TES.tau_el;%[1/s]

% this is the ETF time constant under the assumption of no pole mixing:L->
% infty
Det.TES.tau_etf_simp = Det.TES.tau0 ./(1+ (1-Det.elec.Rl./Det.TES.Ro)./(1+ Det.TES.beta+ Det.elec.Rl/Det.TES.Ro).*Det.TES.LG);%[s]
Det.TES.w_etf_simp   = 1./Det.TES.tau_etf_simp;
%--------------------------------------------------------------------------


% these are the frequencies of the poles taking into account pole mixing:
wp_avg = (1./Det.TES.tau_el)/2+(1./Det.TES.tau_I)/2;
dw= sqrt(((1./Det.TES.tau_el)-(1./Det.TES.tau_I)).^2 - 4.* (Det.TES.Ro ./Det.elec.Lt) .* Det.TES.LG .*(2 + Det.TES.beta)/Det.TES.tau0)/2;%[1/s]

Det.TES.wp_p = wp_avg+dw; % High Frequency Pole ->  ~ L/R
Det.TES.wp_m = wp_avg-dw; % Low Frequency Pole ->  ~ tau_eff


Det.TES.taup_p = 1./Det.TES.wp_p;
Det.TES.taup_m = 1./Det.TES.wp_m;

Det.Response.dIdV0 = (1-Det.TES.LG)./(Det.elec.Rl+Det.TES.Ro.*(1+Det.TES.beta)+ Det.TES.LG.*(Det.TES.Ro-Det.elec.Rl));%[1/Ohm]

%---- Let's check our inversion equations which take wp_p,wp_m,w_z,dVdI(0)
%and get LG,Lt, Tau0, and beta

w_elchk = Det.TES.wp_p+Det.TES.wp_m-Det.TES.w_I;

E= ((w_elchk-Det.TES.w_I).^2-(Det.TES.wp_p - Det.TES.wp_m).^2)./w_elchk./Det.TES.w_I;

beta_chk = 4/(E+4)/Det.Response.dIdV0/Det.TES.Ro-Det.elec.Rl/Det.TES.Ro-1;
LG_chk = (1/Det.Response.dIdV0- (Det.elec.Rl+ Det.TES.Ro.*(1+beta_chk)))./(Det.TES.Ro-Det.elec.Rl+1/Det.Response.dIdV0);
Lt_chk = (Det.elec.Rl+ Det.TES.Ro.*(1+beta_chk))./w_elchk;
tau0_chk= (1-LG_chk)/Det.TES.w_I;

%----- dIdV step function -----
n_t=1e4;
t=linspace(-5,10,n_t)'.*Det.TES.taup_m;%[s]


%let's also caculate the amplitude of the dIdV between the peaks:
dIdVmid = 1./(Det.elec.Rl+Det.TES.Ro.*(1+Det.TES.beta));%[Ohm]

dIdV_chk = Det.Response.dIdV0 .* (1+1j.*omega./Det.TES.w_I)./(1+1j.*omega./Det.TES.wp_p)./(1+1j.*omega./Det.TES.wp_m);

Det.Response.dIdV_step = Det.Response.dIdV0 .*(1 +...
                     -(Det.TES.taup_p-Det.TES.tau_I)./(Det.TES.taup_p-Det.TES.taup_m).*exp(-t./Det.TES.taup_p)+...
                     -(Det.TES.taup_m-Det.TES.tau_I)./(Det.TES.taup_m-Det.TES.taup_p).*exp(-t./Det.TES.taup_m));%[1/Ohm]
Det.Response.dIdV_step(t<0)=0;  
Det.Response.t=t;

if lgc_plt
    %dIdPt Magnitude
    figure(1)
    plot(omega/(2*pi),abs(Det.Response.dIdPt),'-k')
    xlabel('Frequency [Hz]')
    ylabel('dIdP [1/V]')
    title('Magnitude of dIdPt')
    set(gca,'yscale','log','xscale','log')
    grid on
    
    figure(2)
    plot(omega/(2*pi),angle(Det.Response.dIdPt)*180/pi,'-k')
    xlabel('Frequency [Hz]')
    ylabel('Phase of dIdPt [Deg]')
    title('Phase of dIdPt')
    set(gca,'xscale','log')
    grid on
    
    %dIdV Magnitude
    figure(3)
    plot(omega/(2*pi),abs(Det.Response.dIdV),'-k')
    hold on
    %plot(omega/(2*pi),abs(dIdV_chk),'--r')
    plot(omega/(2*pi), dIdVmid .* ones(n_omega,1),'--c')
    plot(omega/(2*pi), abs(Det.Response.dIdV0).* ones(n_omega,1),'--c')
    hold off
    xlabel('Frequency [Hz]')
    ylabel('dIdV [1/\Omega]')
    title('Magnitude of dIdV')
    set(gca,'yscale','log','xscale','log')
    grid on

    figure(4)
    plot(omega/(2*pi),angle(Det.Response.dIdV)*180/pi,'-k')
    hold on
    plot(omega/(2*pi),angle(dIdV_chk)*180/pi,'--r')
    hold off
    xlabel('Frequency [Hz]')
    ylabel('Phase of dIdV [Deg]')
    title('Phase of dIdV')
    set(gca,'xscale','log')
    grid on
     
    % let's look at the dIdV step function
    figure(5)
    plot(t,Det.Response.dIdV_step,'-k')
    hold on
    plot(t, Det.Response.dIdV0 .* ones(n_t,1),'--c')
    plot(t, dIdVmid .* ones(n_t,1),'--c')
    hold off
    xlabel('time [s]')
    ylabel('dIdV Step Function [1/Ohm]')
    title('Step function voltage response')
    grid on
    
    %let's look at the standard dVdI complex vs real plot
    figure(6)
    set(6,'position',[100 100 800 450])
    plot(real(Det.Response.Ztot),imag(Det.Response.Ztot),'-k')
    hold on
    plot(real(Det.Response.Ztes),imag(Det.Response.Ztes),'-b')
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
    
    legend({'Z_{tot}','Z_{tes}'},'location','best')
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