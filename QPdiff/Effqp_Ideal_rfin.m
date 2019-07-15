function [eff_absb,Tk,k]=Effqp_Ideal_rfin(lfin,ltrap,lgc_plt)
%let's calculate the average quasi-particle absorption into the TES for a
% 1D RECTANGLUAR Fin.
%
% All inputs must be the same size
%
% INPUTS:
%           1) Fin Length
%           2) Trapping Length defined as sqrt(D*tau_trap)
%           3) W/Al interface transmission probability
%           4) scattering length in the Al fin
%
% If the interface has perfect transmission then only send the first two
% inputs!

%/////////// Input Check //////////////////////////////////////////////////

if nargin<3
    lgc_plt=false;
end

% First, let's find the input with the largest number of elements:
dim_in = {size(lfin),size(ltrap)};
niter = [numel(lfin),numel(ltrap)];

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

if lgc_bad
    display('INPUTS ARE NOT CORRECTLY FORMATTED!')
    return
end

%Finally let's reshape all of these expanded inputs to a single column vector
if niter > 1
    lgc_reshape=true;
    
    lfin=reshape(lfin,[niter,1]);
    ltrap=reshape(ltrap,[niter,1]);
else
    lgc_reshape=false;
end

%/////////// Find Wave Vectors ////////////////////////////////////////////
%How many basis vectors will we use to approximate the system
nm=60;
m=[0:nm-1];

% The k-values are easy to calculate by hand for perfect transmission
k= (pi./(2*lfin))*(2*m+1);

%/////////// Fit Initial Conditions ///////////////////////////////////////
% The amplitudes for each of the basis functions for a uniform
% distribution are also easily calculated:
Tk= (4*(-1).^m)./(pi*(2*m+1));

if lgc_plt
    %let's plot the initial distribution as a function of x for the
    %first lfin
    nx=1000;

    xfin = linspace(0,lfin(1),nx)';

    Yk = cos(xfin*k(1,:));

    TYk = Yk .* (ones(nx,1)*Tk);
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
%//////////// Calculate Collection Efficiency /////////////////////////////
mLG= ones(niter,1)*m;
ltrapLG= ltrap*ones(1,nm);

eff_absb= sum(8./(pi.^2.*(2.*mLG+1).^2) .* (1- 1./((ltrapLG.*k).^2+1)),2);

%//////////// Reshape Ouput ///////////////////////////////////////////////
%if we reshaped the input arrays ... let's return the correctly dimensioned eff_absb
if lgc_reshape
    eff_absb=reshape(eff_absb,dim_in);
end    

