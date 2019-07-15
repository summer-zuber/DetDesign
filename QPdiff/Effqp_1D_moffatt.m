function [fQP,ld,la,eff_absb] = Effqp_1D_moffatt(lfin,hfin,loverlap,eff_absb)
% Here we're using Robert Moffatt's full QP model. There are some pretty
% big assumptions:
%   1) diffusion length scales with Al thickness (trapping seems to be
%   surface dominated! and diffusion seems to be thickness limited)
%   2) 
%
%  Future:  scale boundary impedance with W thickness
%
%  INPUTS:
%   1) fin length [um]
%   2) fin height [um]
%   3) W/Al overlap [um]
%   4)W/Al transmission/trapping probability (Default = 1.22e-4)
%
%  OUTPUTS:
%   1) quasi-particle collection efficiency
%   2) diffusion length for Al
%   3) W/Al surface absorption length
%   4) W/Al transmission probability
%
%
% From Robert Moffatt's pdf:  Analytical_banana_v3.pdf (attached to Robert Moffatt's email pdf:  15/06/07 "talk today or tomorrow")
%                             https://www.stanford.edu/~rmoffatt/Papers/Diffusion%20and%20Absorption%20with%201D%20and%202D%20Solutions.pdf
%                             DOI: 10.1007/s10909-015-1406-7
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    % fractional trapping/transmission coefficient at W/Al interface
    %eff_absb = 1e-4; %->  this matches rmoffat's email
    eff_absb = 1.22e-4; % -> this matches J. Yen et al paper 
    %eff_absb=1e-2;
end

%vqp=6e5; %[m/s]
%Dqp = vqp*hfin; %[m^2/s]

%---- pertinent length scales (fits from Stanford QP data) ----
% diffusion length
%ld = 666.*hfin;%[um]->  this was Jeff's original fit
ld = 567.*hfin;%[um] -> this is the fit in Jeff's published LTD 16 data
% surface impedance "length"
la = 1/eff_absb*hfin.^2 ./ loverlap; %[um]

%--- dimensionless quantities ---
% dimensionless diffusion length
LambdaD = ld./lfin;
% dimensionless surface impedance
lambdaA = la./lfin;

%let's convert lambdaA into a intuitive dimensionless quantity
% this is basically the probability for transmission
eff_absb = hfin.^2 ./(loverlap.*la);


% Quasi-Particle collection coefficient
fQP = LambdaD.^2 ./ (LambdaD.* coth(1./LambdaD)+lambdaA);