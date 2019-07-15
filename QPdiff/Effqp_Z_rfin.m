function eff_absb=Effqp_Z_rfin(lfin,ltrap,fwal,lscat,lgc_plt)
%let's calculate the average quasi-particle absorption into the TES for a
% 1D RECTANGLUAR Fin with an impedance boundary condition.
%
% All inputs must be either 1x1 or the same size ndimensional array
%
% INPUTS:
%           1) Fin Length
%           2) Trapping Length defined as sqrt(D*tau_trap)
%           3) W/Al interface transmission probability
%           4) scattering length in the Al fin
%           5) lgc_plt
%
% If the interface has perfect transmission then only send the first two
% inputs!

%How many basis vectors will we use to approximate the system
nm=60;
m=[0:nm-1];

%/////////// INPUT QUALITY CHECKS /////////////////////////////////////////

if nargin < 5
    lgc_plt=false;
end

% First, let's find the input with the largest number of elements:
dim_in = {size(lfin),size(ltrap),size(fwal),size(lscat)};
niter = [numel(lfin),numel(ltrap),numel(fwal),numel(lscat)];

[niter,ind] = max(niter);
dim_in = dim_in{ind};

%Now let's go through and check to make sure that  all variables either
%have this dimension or are [1,1]

lgc_bad=false;

%lfin
if  prod(size(lfin))==1
    lfin=repmat(lfin,dim_in);
elseif any(~(size(lfin)==dim_in))
    lgc_bad=true;
end

%ltrap
if  prod(size(ltrap))==1
    ltrap=repmat(ltrap,dim_in);
elseif any(~(size(ltrap)==dim_in))
    lgc_bad=true;
end

%fwal
if  prod(size(fwal))==1
    fwal=repmat(fwal,dim_in);
elseif any(~(size(fwal)==dim_in))
    lgc_bad=true;
end

%lscat
if  prod(size(lscat))==1
    lscat=repmat(lscat,dim_in);
elseif any(~(size(lscat)==dim_in))
    lgc_bad=true;
end

if lgc_bad
    display('INPUTS ARE NOT CORRECTLY FORMATTED!')
    return
end

%Finally let's reshape all of these expanded inputs to a single column vector
if niter > 1
    lgc_reshape=true;
    
    lfin=reshape(lfin,[niter,1]);
    ltrap=reshape(ltrap,[niter,1]);
    fwal=reshape(fwal,[niter,1]);
    lscat=reshape(lscat,[niter,1]);
else
    lgc_reshape=false;
end

%let's initialize output 
eff_absb=ones(niter,1);

for jiter=1:niter
    %/////////////// Find Wave Vectors (k) ////////////////////////////////////

    % The boundary conditions on the two sides quantize the allowed k vectors.
    % Unlike in the ideal transmission case, the quantization equation must be
    % solved computationally

    %FYI:kwal in my notes is called beta
    kwal= 3/2.*fwal(jiter)./lscat(jiter);

    k=ones(nm,1);
    for jm=1:nm
        %let's guess a k-value for the zero!
        ko = pi./(2*lfin(jiter)).*2*m(jm);
        kguess = abs(real(ko + 1i*log((1i*ko+kwal)/(1i*ko-kwal))/(2*lfin(jiter))));

        k(jm) = fzero(@kdef_Zboundary,kguess,[],lfin(jiter),kwal,m(jm));
    end    

    %////////////// Match Initial Conditions //////////////////////////////////
    % To calculate the amplitudes for the basis functions we must least
    % sqaures fit the problem!

    %let's fit the system at 3*nm evenly spaced locations along the fin
    nfit=3*nm;
    xfin= linspace(0,lfin(jiter),nfit)';

    %We're solving pde for the case that the particles are smoothly
    %distributed
    n0=ones(nfit,1);

    %spatially the different k modes look like cosines
    Yk= cos(xfin* k');

    Tk=lsqlin(Yk,n0);

    if lgc_plt & mod(jiter,10)==1
        %let's plot the initial distribution as a function of x for the
        %first lfin
        nx=1000;

        xfin = linspace(0,lfin(jiter),nx)';

        Yk = cos(xfin*k');

        TYk = Yk .* (ones(nx,1)* Tk');
        Ytot = sum(TYk,2);

        cmap=colormap(jet(nm));

        figure(1)
        clf(1)
        hold on
        for jm=1:nm
            plot(xfin/lfin(1),TYk(:,jm),'-','color',cmap(jm,:))
        end
        plot(xfin/lfin(1),Ytot,'-k')
        hold off

        xlabel('x/l_{fin}')
        ylabel('n(x)')
        title('Initial Density Distribution')
        drawnow
    end

    % The trapping efficiency can now be calculated  
    eff_trap = sum(Tk .* sin(k.*lfin(jiter))./(k.*lfin(jiter)) ./ ((ltrap(jiter).*k).^2+1));

    % The collection efficiency is simply 1-trapping efficiency
    eff_absb(jiter) = 1-eff_trap;
end %jiter

%if we reshaped the input arrays ... let's return the correctly dimensioned eff_absb
if lgc_reshape
    eff_absb=reshape(eff_absb,dim_in);
end    

end %function

function delk = kdef_Zboundary(k,lfin,kwal,m);
%this is the function which defines the quantized k values for the surface
%impedance
    delk =k- pi./(2*lfin).*(2*m + 1i*log((1i*k+kwal)/(1i*k-kwal))/pi);
    delk=real(delk);
end
