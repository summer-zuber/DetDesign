function [R_parallel] = Rtes_parallel(R)
% This function calculates the total resistance of a bunch of TES 

R_parallel = 1./sum(1./R);

%Det.TES.Rtes_parallel = 1./sum(1./Det.TES.R);



end

