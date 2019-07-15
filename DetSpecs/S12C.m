function [Det]= S12C()
% Properties of S12C, an iZIP4 detector
%
% - We're going to treat each side as a Separate Detector (just multiply the channel resolution by sqrt(8))
% 
%
% Created: 12/3/13 MCP
%--------------------------------------

% let's inherit all the iZIP4 specs;
[Det]=iZIP4(fUCB(),eCDMSII(fUCB),'Si','W','Al');


% now let's replicate 8 times to keep track of all the channels
Det=repmat(Det,[2,1]);

%//////////////// Measured Quantities /////////////////////
NameStr={'S12C low Tc'
         'S12C High Tc'};

% Tc
Tc=[65.5,100.5]*1e-3;%[K]
 
% Baseline energy resolution
% No low frequency noise subtraction
resPt_raw = [102,175]/2.35;%eV;%[sigma eVt]
% low frequency noise subtraction
resPt_bs = [55,143]/2.35;%eV;%[sigma eVt]

% Fraction of phonon energy absorbed and dissappated by TES
fPtAbsb = 0.12 ; %the systematics are significant on this fraction usually

for jj=1:2
    Det(jj).name=NameStr{jj};
    Det(jj).Tc=Tc(jj);
    Det(jj).resPt_raw=resPt_raw(jj);
    Det(jj).resPt_bs=resPt_bs(jj);
    Det(jj).fPtAbsb=fPtAbsb;
end    