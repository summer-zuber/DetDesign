function [Det]= G3D();
% Properties of G3D, an iZIP2 detector
%
% Created: 12/3/13 MCP
%--------------------------------------

% let's inherit all the iZIP2 specs
Det= iZIP2(fUCB(),eCDMSII(),'Ge');

%//////////////// Measured Quantities /////////////////////
% name
Det.name = 'G3D';
% Tc
Det.TES.Tc=65e-3;%[K]

% Baseline energy resolution
Det.resPt = 32;%[sigma eVt]

% Fraction of phonon energy absorbed and dissappated by TES
Det.fPtAbsb = 0.2 ; %the systematics are significant on this fraction usually