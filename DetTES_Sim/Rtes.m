function [R,alpha,beta] = Rtes(I,T,Det)
% This function finds the resistance of a TES from the current and the temperature.
%
% Note: This function takes both vectors and scalars ... depending upon the
% calling function, this one function will handle (TES is parallel, TES in
% series, or a single simple TES)
%
% 17/03/30: Created MCP
% -------------------------------------------------------------------------

%cross-sectional area for each "simulated" TES (in the simulation we group physical TES into groups to make the matrices tractable)
Acs = (Det.TES.n * Det.TES.t *Det.TES.w);
Ai =Det.TES.dAidAcs .* Acs;
    

chi= (T-Det.TES.Tc +  (I/Ai).^(2/3) )./Det.TES.wTc;


R = Det.TES.Rn/2 .*(1.0 + tanh(chi));

dRdchi=Det.TES.Rn/2 .*(1.0 - tanh(chi).^2);
dRdT= dRdchi ./Det.TES.wTc;
dRdI= 2/3*dRdchi./(I.^(1/3) .* Det.TES.wTc .* Ai.^(2/3));

alpha = dRdT.*T./R;
beta  = dRdI.*I./R;
