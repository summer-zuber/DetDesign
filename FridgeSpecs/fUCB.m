function [fUCB] = fUCB(T_MC)
if nargin<1 | isempty(T_MC)
    T_MC=35e-3;
end
    
% Fridge specification object which corresponds to 75uW at UCB 
fUCB=[];
fUCB.name='UCB';
fUCB.T_MC = T_MC; %[K]
fUCB.T_Still = 1.4; %[K]
fUCB.T_4K = 4.0; %[K]


end

