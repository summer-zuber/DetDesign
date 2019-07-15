function [eff_absb]=Effqp_absb_rfin(lfin,ltrap);
%let's calculate the average quasi-particle absorption assuming perfect
%collection at the interface for a 1-Dimensional rectangular fin.

nj=40; % highest k-mode

%let's assume that lfin and ltrap are both large arrays
eff_absb=0;
for m=2*[0:nj]+1
    eff_absb = eff_absb + (ltrap./lfin).^2 ./(1+(ltrap./lfin).^2.*m.^2 .* pi.^2 /8);
end

%let's have a dummy code as well ... just for kicks
%eff_absb= exp(-lfin./ltrap);
