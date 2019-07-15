function Det = Equilibrium_1tes(Det)
%
% To find the equilbrium point we solve the set of non-linear equations in dITVb_Equilibrium_1TES:
%
% 1)  L dIdt = Vb - I (Rl +R)
% 2)  C dTdt = I^2 R - Pcool
% 3)     0   = R - Rwanted
%
% X = Equilibrium Point Phase Space Vector
%     I ->  TES current 
%     T ->  TES temperature
%     Vb ->  Bias Voltage
%
%--------------------------------------------------------------------------

%Initial Guess:
% We find the equilibrium point assuming that beta = 0;

Det=SimpleEquilibrium_1TES(Det);

Xguess = [Det.TES.Io;Det.TES.To;Det.TES.Vbias];

%let's hide all the other dependencies in function call
dY  = @(X) dITVb_Equilibrium_1TES(X,Det);
Xeq = fsolve(dY,Xguess);

Det.TES.Io=Xeq(1);
Det.TES.To=Xeq(2);
Det.TES.Vbias=Xeq(3);


%\\\\\\\\\\\ Equilibrium Point Properties \\\\\\\\\\\\\
pc=PhysicalConstants();

%------- Alpha/Beta/Loop Gain --------------------------------------------

[Det.TES.Ro,Det.TES.alpha,Det.TES.beta]=Rtes(Det.TES.Io,Det.TES.To,Det);

Det.TES.LG= Det.TES.alpha/Det.TES.nPep*(1- (Det.fridge.T_MC/Det.TES.To).^Det.TES.nPep);
  
%------ TES Properties @ Equilibrium --------------------------------------------------
Det.TES.C   =  (Det.TES.fCsn*Det.TES.gC_v.* Det.TES.To)*(Det.TES.vol); %J/K
Det.TES.Gep =  Det.TES.nPep*Det.TES.gPep_v*(Det.TES.vol).*  Det.TES.To.^(Det.TES.nPep-1); %W/K
Det.TES.Po  =               Det.TES.gPep_v*(Det.TES.vol).* (Det.TES.To.^(Det.TES.nPep)-Det.fridge.T_MC.^(Det.TES.nPep)); %W
Det.TES.Io  = sqrt(Det.TES.Po/Det.TES.Ro);%[A] 
Det.TES.tau0= Det.TES.C./Det.TES.Gep; %[s]
