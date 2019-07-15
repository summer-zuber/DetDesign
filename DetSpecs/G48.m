function [Det]= G48()
% Properties of G48, an iZIP4 detector
%
% - We're going to treat each channel as a Separate Detector (just multiply the channel resolution by sqrt(8))
% 
%
% Created: 12/3/13 MCP
%--------------------------------------

% let's inherit all the iZIP4 specs
[Det]=iZIP4(fUCB(),[],'Ge','W','Al');


% now let's replicate 8 times to keep track of all the channels
Det=repmat(Det,[8,1]);

%//////////////// Measured Quantities /////////////////////
NameStr={'G48 As1'
         'G48 Bs1'
         'G48 Cs1'
         'G48 Ds1'
         'G48 As2'
         'G48 Bs2'
         'G48 Cs2'
         'G48 Ds2'};

% Tc
Tc=[88 91 91 91 109 101 101 101]*1e-3;%[K]
 
% Baseline energy resolution
% http://titus.stanford.edu/cdms_restricted/detector_physics/iZIP/ebook/101008/Pnoise_g48.html
resPt_bs = [34 35 38 33 39 38 37 46]*sqrt(8);%[sigma eVt]
resPt_raw = [66 66 66 240 60 78 43 132]*sqrt(8); %[sigma eVt]
% Fraction of phonon energy absorbed and dissappated by TES
fPtAbsb = 0.12 ; %the systematics are significant on this fraction usually

for jj=1:8
    Det(jj).name=NameStr{jj};
    Det(jj).Tc=Tc(jj);
    Det(jj).resPt_raw = resPt_raw(jj);
    Det(jj).resPt_bs = resPt_bs(jj);
    Det(jj).fPtAbsb = fPtAbsb;
end    