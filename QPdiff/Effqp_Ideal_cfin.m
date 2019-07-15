function [eff_absb]=Effqp_absb_cfin(lfin,lwal,ltrap,lgc_plt)
%let's calculate the average quasi-particle absorption into the TES for a
% 2D cylindrical fin WITH PERFECT ABSORPTION ON THE TES INTERFACE:
%
% All inputs must be either 1x1 or a single arbitrary m-dimensional array
%
% INPUTS:
%           1) Fin Length
%           2) W/A total overlap length (this is usually between ltes & 2ltes)
%           3) Trapping Length defined as sqrt(D*tau_trap)
%           4) W/Al interface transmission probability
%           5) scattering length in the Al fin
%

%//////////// INPUT CHECKS ////////////////////////////////////////////////

if nargin==0
    lfin = 2.0;
    lwal = 0.5;
    ltrap= 1.0;
end    

if nargin<4
    lgc_plt=false;
end

% First, let's find the input with the largest number of elements:
dim_in = {size(lfin) , size(lwal) ,size(ltrap) };
niter  = [numel(lfin), numel(lwal),numel(ltrap)];

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

%lwal
if  prod(size(lwal))==1
    lwal=repmat(lwal,dim_in);
elseif any(~(size(lwal)==dim_in))
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
    
    lfin = reshape(lfin,[niter,1]);
    lwal = reshape(lwal,[niter,1]);
    ltrap= reshape(ltrap,[niter,1]);
else
    lgc_reshape=false;
end

%let's initialize output 
eff_absb=ones(niter,1);

for jiter=1:niter
    if mod(jiter,10)==1
        display(['Effqp_Ideal_cfin: event number ',num2str(jiter),' of ',num2str(niter)])
    end    
    %How many basis vectors will we use to approximate the system
    nk=80;

    %The collesium QET design isn't actually a perfectly azimuthally symmetric
    %system ... this is definitely an abstraction. So for a decent guess of the
    %inner radius, let's just say that the total amount of overlap is equal to
    %the circumference of a circle with radius Rin.

    Rin  = lwal(jiter)/(2*pi);

    %The outer radius of they system is even easier to approximate ... it's
    %just Rin+ the fin length
    Rout = Rin + lfin(jiter);

    %//////////// k vector quantization ///////////////////////////////////////

    % We need to find the wave vectors which satisfy the boundary conditions:
    %   1) @Rin -> perfect transmission
    %   2) @Rout -> perfect reflection

    % Just by plotting the quantizing condition in matlab we've found that the
    % k's are very, very similar to those of the cartesian system! So let's use
    % the cartesian system as the first guess.
    k0guess= pi/(Rout-Rin)*(2*[0:nk-1]'+1)/2; 

    k0_low  = pi/(Rout-Rin)*[0:nk-1]';
    k0_high = pi/(Rout-Rin)*[1:nk]';

    %let's modify the lower bound to be 5% of the upper bound for the 0 basis
    %vector... for the Neumann functions z=0 is infinite!
    k0_low(1)= 0.05*k0_high(1);


    k0= zeros(nk,1);
    for jk=1:nk
        %k0(jk)=fzero(@matchBC,k0guess(jk),[],Rin,Rout);
        k0(jk)=fzero(@matchBC,[k0_low(jk),k0_high(jk)],[],Rin,Rout);
    end  

    if lgc_plt & mod(jiter,10)==1

        %let's plot k0guess, k0, and the quantizing condition ... we can use
        %this to make certain that fzero didn't miss any zeros.
        
        dk1D= pi/(Rout-Rin);
        
        k= linspace(0,k0guess(end),1e5)';

        delq= matchBC(k,Rin,Rout);
        
        h=zeros(2,1);
        figure(1)
        plot(k/dk1D,delq,'k')
        hold on
        h(1)=plot(k0guess/dk1D,zeros(nk,1),'*b')
        h(2)=plot(k0/dk1D,zeros(nk,1),'*g')
        hold off
        xlabel('k / (\pi/(r_{o}-r_{i})) ')
        title('k Quantization Constraint for Cylindrical Boundary Conditions')
        legend(h,{'1D k guesses','2D k solutions'})
        xlim([0,15])
        ylim([-.4,.4])
    end

    %///////////// Calculate Jo/Yo Mixing  ////////////////////////////////////

    %next let's find the mixing vector for the J and Y bessel components:
    Fjy= [ bessely(1,k0*Rout), -besselj(1,k0*Rout)];

    %let's normalize Fjy (just because):
    Fjy= Fjy./ sqrt(sum(Fjy.^2,2)*ones(1,2));

    if lgc_plt
       %let's plot all the basis functions in x-space just to prove that everyting is working


        nr=1e4;
        r= linspace(Rin,Rout,nr)';

        Yk= zeros(nr,nk);
        for jk=1:nk
            Yk(:,jk)= Fjy(jk,1)*besselj(0,k0(jk)*r)+ Fjy(jk,2)*bessely(0,k0(jk)*r);
        end

        cmap=colormap(jet(nk));
        figure(2)
        clf(2)
        hold on
        for jk=1:nk
            plot(r,Yk(:,jk),'color',cmap(jk,:))
        end
        hold off
        xlabel('R')
        title('Basis Functions for QP diffusion in cylindrical fins')
    end

    %\\\\\\\\\\\\\\\\\\ Initial Condition Matching \\\\\\\\\\\\\\\\\\\\\\\\\\\\

    % Now that we know the basis function which solve this problem, let's least
    % squares fit to the initial distribution to find the coefficients

    %let's over constrain the fit by look at 3*nk evenly spaced locations along the fin
    nfit=3*nk;

    r= linspace(Rin,Rout,nfit)';

    Yk= zeros(nfit,nk);
    for jk=1:nk
        Yk(:,jk)= Fjy(jk,1)*besselj(0,k0(jk)*r)+ Fjy(jk,2)*bessely(0,k0(jk)*r);
    end

    %the initial density distribution is totally even:
    n0 = ones(nfit,1);

    Tk = lsqlin(Yk,n0);

    if lgc_plt
         % Let's plot the initial density distribution:

        nr=1e4;
        r= linspace(Rin,Rout,nr)';

        Yk= zeros(nr,nk);
        for jk=1:nk
            Yk(:,jk)= Fjy(jk,1)*besselj(0,k0(jk)*r)+ Fjy(jk,2)*bessely(0,k0(jk)*r);
        end

        Ytot= Yk*Tk;

        cmap=colormap(jet(nk));
        figure(3)
        clf(3)
        hold on
        for jk=1:nk
            plot(r,Tk(jk)*Yk(:,jk),'--','color',cmap(jk,:))
        end
        plot(r,Ytot,'-k')
        hold off
        xlabel('R')
        title('Homogenous QP Density Initial Condition Fit')

    end

    %///////// Trapping Efficiency ////////////////////////////////////////////

    %let's calculate the total number of particles at t=0

    N0 =        sum((2*pi* Tk)./(k0.^2).*(Fjy(:,1).*(k0.*Rout.*besselj(1,k0.*Rout)-k0.*Rin.*besselj(1,k0.*Rin))...
                                         +Fjy(:,2).*(k0.*Rout.*bessely(1,k0.*Rout)-k0.*Rin.*bessely(1,k0.*Rin))...
                                         ));

    eff_absb(jiter) =  sum((2*pi* Tk)./(k0.^2).*(ltrap(jiter).^2 .* k0.^2)./(ltrap(jiter).^2 .* k0.^2+1)...
                    .*( Fjy(:,1).*(k0.*Rout.*besselj(1,k0.*Rout)-k0.*Rin.*besselj(1,k0.*Rin))...
                       +Fjy(:,2).*(k0.*Rout.*bessely(1,k0.*Rout)-k0.*Rin.*bessely(1,k0.*Rin))...
                      ))./N0;
end % for jiter

%//////////// Reshape Ouput ///////////////////////////////////////////////
%if we reshaped the input arrays ... let's return the correctly dimensioned eff_absb
if lgc_reshape
    eff_absb=reshape(eff_absb,dim_in);
end

end % function

function delq = matchBC(k,ri,ro);
%this is the function which defines the quantized k values for a fin with
%azimuthal symmetry. The boundary conditions are:
%
%   1) at ri -> perfectly absorbing boundary
%   2) at ro -> perfectly reflecting boundary

    delq= bessely(0,k*ri).*besselj(1,k*ro)-besselj(0,k*ri).*bessely(1,k*ro);

    
end
