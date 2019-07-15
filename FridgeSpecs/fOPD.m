function [fOPD] = fOPD()
% Fridge specification object which corresponds nominally to that needed
% for the ultimate optical phonno detector experiment
fOPD=[];
fOPD.name='OPD';
fOPD.T_MC = 7e-3; %[K] from 2013 G2 SNOLAB Proposal
fOPD.T_CP = 50e-3; %[K] from 2013 G2 SNOLAB Proposal
fOPD.T_Still=900e-3;%[K] from 2013 G2 SNOLAB Proposal
fOPD.T_4K = 3.8; %[K] from 2013 G2 SNOLAB Proposal
fOPD.Qparasitic=0;%[W]
end

