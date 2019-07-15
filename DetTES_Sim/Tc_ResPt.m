function [ResPt]=Tc_ResPt(Tc,DetFunc,lgcPlt,fridge,elec,MatAbsb)
% In this script we calculate the estimated energy sensitivity of a
% detector as a function of Tc
%
% inputs:
%   1) Tc -> a vector [K]
%   2) DetFunction -> a function handle: example  DetFunc  = @iZIP7
%                   Default = @iZIP7
%   3) lgcPlt -> true/false 
%               (Default=false)
%   FOR MORE INFORMATION LOOK UP "function handle" in matlab help
%   4) fridge -> fridge object
%                default: [] -> goes to the default for that detector type
%   5) elec ->  electrics specs
%                default: [] -> goes to the default
%   6) Absorber Material 
%                default: []-> goes to the default for that detector type
%                'Si'
%                'Ge'
% 12/8/13 MCP
%--------------------------------------------------------------------------

if nargin<1 | isempty(Tc)
    Tc= logspace(log10(10e-3),log10(100e-3),30)';
end    
nTc=length(Tc);

if nargin<2 | isempty(DetFunc)
    DetFunc = @iZIP7;
end
%let's check to make certain that DetFunc is a function handle
if ~isa(DetFunc,'function_handle')
    display('Need to use a function handle in input slot 2!')
    return
end

if nargin<3
    lgcPlt=false;
end
 
if nargin<4 | isempty(fridge)
    %let's use the default chosen in the specific device
    fridge=[];
end  

if nargin<5 | isempty(elec)
    %let's use the default chosen in the specific device
    elec=[];
end  


if nargin<6 | isempty(MatAbsb)
    %let's use the default chosen in the specific device
    MatAbsb=[];
end    

MatTES = MaterialProperties('W');

ResPt=zeros(nTc,1);
for jj=1:nTc
    % let's change the Tc parameter in the W object:
    MatTES.Tc = Tc(jj);
    
    %let's scale the transition width with Tc
    wTc_1090 =1.4e-3* MatTES.Tc/68e-3; %[K]
    MatTES.wTc = wTc_1090/2/log(3); %[K]
    
    % Create the detector object:
    %HV4mm_4ch(fridge,elec,MatAbsb,MatTES,MatQET,ltes,lfin,hfin,loverlap,Tc,wTc,T_MC,gPep_v)
    Det = DetFunc(fridge,MatAbsb,MatTES);
    
    if Det.TES.Tc > Det.fridge.T_MC*1.01 
        % Now let's calculate the energy resolution:
        ResPt(jj) = SimulatedNoise_1TES(Det,lgcPlt);
    else
        ResPt(jj) = NaN;
    end    
end


%---------------- plot -------------------------------

figure(1)
plot(Tc*1e3,ResPt,'-k')
%hold on

xlim([min(Tc),max(Tc)]*1e3)
ylim([0.9*min(ResPt),1.1*max(ResPt)])
grid on
set(gca,'yscale','log','xscale','log')

title([Det.name,': Energy Resolution vs TES Tc'])
ylabel('Phonon Energy Resolution (sigma)  [eVt]')
xlabel('TES Tc [mK]')



