function [fQP,ld,la,eff_absb] = Effqp_2D_moffatt(lfin,hfin,loverlap,ltes,eff_absb);
% Here we're using Robert Moffatt's full QP model. There are some pretty
% big assumptions:
%   1) diffusion length scales with Al thickness (trapping seems to be
%   surface dominated! and diffusion seems to be thickness limited)
%   2) 
%
%  Future:  scale boundary impedance with W thickness
%
%  INPUTS:
%   1) fin length   [um]
%   2) fin height   [um]
%   3) W/Al overlap [um]
%   4) TES length   [um]
%   5) W/Al transmission/trapping probability: (Default 1.22e-4)
%
%  OUTPUTS:
%   1) quasi-particle collection efficiency
%   2) diffusion length
%   3) W/Al surface absorption length
%   4) W/Al transmission probability
%
% From Robert Moffatt's email pdf:  15/06/07 "talk today or tomorrow"
%                                   https://www.stanford.edu/~rmoffatt/Papers/Diffusion%20and%20Absorption%20with%201D%20and%202D%20Solutions.pdf
%                                   DOI: 10.1007/s10909-015-1406-7
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    % fractional trapping/transmission coefficient at W/Al interface
    %eff_absb = 1e-4; %->  this matches rmoffat's email
    eff_absb = 1.22e-4; % -> this matches J. Yen et al paper 
    %eff_absb=1e-2;
end

%vqp=6e5; %[m/s]
%Dqp = vqp*hfin; %[m^2/s]
%let's first calculate the inner and outer circumference and effective radius for our pie shaped QP
%collection fins
ci = 2.*ltes;
ri = ci./(2*pi);

%let's calculate the outer circumference of a "very simplified" ellipse 
co1 = 2*ltes+2*pi*lfin;
ro1 = co1/(2*pi);

%let's do another approximation
a = lfin+ltes/2;
b = lfin;

% https://www.mathsisfun.com/geometry/ellipse-perimeter.html
h= (a-b).^2./(a+b).^2;
co = pi*(a+b).*(1+3*h./(10+sqrt(4-3*h)));
ro = co1/(2*pi);

%--- pertinent length scales ----
%diffusion length
%ld = 666.*hfin;%[um]->  this was Jeff's original fit
ld = 567.*hfin;%[um] -> this is the fit in Jeff's published LTD 16 data
%surface impedance length
la = 1/eff_absb*hfin.^2./loverlap; %[um] -> these match those values used by Noah
la_chk = (1e6*1600/900^2*5)*hfin.^2 ./loverlap; %[um] 

%---- dimensionless scales ----
rhoi = ri./ld;
rhoo = ro./ld;
lambdaA = la./ld;

% Quasi-Particle collection coefficient
fQP = 2.*rhoi./(rhoo.^2-rhoi.^2)...
      .* (besseli(1,rhoo).*besselk(1,rhoi)-besseli(1,rhoi).*besselk(1,rhoo))...
      ./ ( besseli(1,rhoo).*(besselk(0,rhoi)+lambdaA.*besselk(1,rhoi))+...
          (besseli(0,rhoi)-lambdaA.*besseli(1,rhoi)).*besselk(1,rhoo) );