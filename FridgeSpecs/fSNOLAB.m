function [fSNOLAB] = fSNOLAB()
% Fridge specification object which corresponds nominally to the proposed
% SNOLAB fridge
fSNOLAB=[];
fSNOLAB.name='SNOLAB';
fSNOLAB.T_MC = 20e-3; %[K] from 2013 G2 SNOLAB Proposal
%fSNOLAB.T_MC = 7e-3; %[K] from 2013 G2 SNOLAB Proposal
fSNOLAB.T_CP = 145e-3; %[K] from 2013 G2 SNOLAB Proposal
fSNOLAB.T_Still=900e-3;%[K] from 2013 G2 SNOLAB Proposal
fSNOLAB.T_4K = 4.8; %[K] from 2013 G2 SNOLAB Proposal
fSNOLAB.Qparasitic=0;%[W]
end

