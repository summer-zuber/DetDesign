function [EstDet] = SimplePtScale_Exp(Test, Design);
%let's scale experimental measurements from a test device to that for the
%design.

EstDet=Design;

EstDet.Meas.sigPt0= Test.Meas.sigPt0 .*(Test.Meas.fPt/Design.Meas.fPt) ...
                               .*sqrt((Design.TES.vol.*Design.nP)/(Test.TES.vol.*Test.nP))....
                               .*sqrt(Design.Meas.tau_Pabsb/Test.Meas.tau_Pabsb);