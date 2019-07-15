function F=Ftfn(Tl,Th,n,lgcB)
%this function estimates the noise suppression in a thermal conductance
%because the two sides aren't the same temperature
% inputs:
%       Tl = low temperature side of thermal conductance
%       Th = high temperature side of thermal conductance
%       n  = heating power law for thermal conductance: Q = K(Th^n-Tl^n)
%       lgcB:  true  = ballistic
%              false = diffusive
%
%---------------------------------

%let's expand out Tl and Th so that they are the same size
nl = length(Tl);
nh = length(Th);

if nl == nh
    %everything is good
elseif nl==1 && nl < nh
    %let's expand out Tl
    Tl = Tl*ones(size(Th));
elseif nh==1 && nh < nl
    %let's expand out Th
    Th = Th*ones(size(Tl));
else
    display('something wrong with Ftfn')
    return 
end    
%let's check to make sure that Tl and Th haven't been switched:
lgc = Tl > Th;
if any(lgc)
    Tx = Th;    
    Th(lgc) = Tl(lgc);
    Tl(lgc) = Tx(lgc);
end

if lgcB
    F = ((Tl./Th).^(n+1)+1)./2;
else    
    F = n/(2*n+1).* ((Tl./Th).^(2*n+1)-1)./((Tl./Th).^(n)-1);
end
end