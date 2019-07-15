% Here we estimate the efficiency on phonon energy collection as a function
% of QET fin geometries

function [QET,B1997]= effQP_simple(lfin,h,lwal)

if nargin < 1
    lfin = 250e-6; %[m]
end
if nargin < 2
    h = 300e-9;%[m]
end
if nargin <3
    lwal =3e-6; %[m]
end

%average QP velocity
vqp = 6e5; %[m/s] -> From Blas' pdf (should probably check)

%------- let's scale everything from the original 1997 banana ------------
B1997=[];
% fin length
%we divide by x2 to turn these numbers into those that match a standard QET fin (TES on only one side)
B1997.lfin= 350e-6/2;%[m]

%banana Al film thickness
B1997.h = 150e-9; %[m]

%banana overlap
B1997.lwal = 3e-6; %[m]

%banana diffusion length
B1997.ltrap = 180e-6;%[m]

%surface transport ->  assuming 1D geometry
%B1997.fs = 1/400; %[1]
%surface transport -> switching to actual geometry
B1997.fs =1/400*(B1997.h/B1997.lwal);

%let's assume thickness limited transport:
B1997.D = B1997.h * vqp/3; %[m^2/s]

%let's calculate tau_trap
B1997.tautrap = B1997.ltrap.^2/B1997.D;%[s]

% Now let's make the various resistances [s/m^2]
%QP trapping resistance
B1997.Rtrap = B1997.tautrap / (B1997.lfin*B1997.h); %[s/m^2]

%QP surface resistance
B1997.Rs  = 1/(B1997.fs*vqp*lwal);%[s/m^2]

%QP diffusive impedance
B1997.Rd  = B1997.lfin/(B1997.D*B1997.h);%[s/m^2]

B1997.Eff = B1997.Rtrap/(B1997.Rtrap+B1997.Rs+B1997.Rs);
%--------------------------------------------------------------------------
%let's scale from B1997
scale = B1997;

QET=[];

QET.lfin = lfin;
QET.h = h;
QET.lwal = lwal;

%let's asssume that D is thickness limited!
QET.D = QET.h * vqp/3; %[m^2/s]

%let's assume QP trapping on the surface (scales with thickness of Al film)
QET.tautrap = scale.tautrap* QET.h/scale.h; %[s]

%let's assume that surface impedance scales linearly with the film
%thickness
QET.fs = scale.fs * QET.lwal/scale.lwal;

QET.ltrap = sqrt(QET.D*QET.tautrap);%[m]

% Now let's make the various resistances [s/m^2]
%QP trapping resistance
QET.Rtrap = QET.tautrap / (QET.lfin*QET.h); %[s/m^2]

%QP surface resistance
QET.Rs  = 1/(QET.fs*vqp*lwal);%[s/m^2]

%QP diffusive impedance
QET.Rd  = QET.lfin/(QET.D*QET.h);%[s/m^2]

QET.Eff = QET.Rtrap/(QET.Rtrap+QET.Rs+QET.Rs);

