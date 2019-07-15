function [Q_TES,Q_Rs,Q_Rp,Tc] = Tc_Heat(Tc,DetFunc,lgcPlt,fridge,MatAbsb,MatTES,MatQET)
% In this script we calculate the fridge heating loads from the phonon
% channels of a detector
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
%   5) MatAbsb -> Object that holds material Properties of Absorber
%                default: []
%   6) MatTES -> Object that holds material properties of TES
%                default: W
%   7) MatQET -> Object that holds material properties of QET
%                default: Al
%
% 12/8/13 MCP
%--------------------------------------------------------------------------

if nargin<1 | isempty(Tc)
    %Tc= logspace(log10(10e-3),log10(100e-3),50)';
    Tc= [20:100]'*1e-3;%[K]
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

if nargin < 3 | isempty(lgcPlt)
    lgcPlt=false;
end    
 
if nargin<4 | isempty(fridge)
    %let's use the default chosen in the specific device
    fridge=[];
end  

if nargin<5 | isempty(MatAbsb)
    %let's use the default chosen in the specific device
    MatAbsb=[];
end    

if nargin < 6 | isempty(MatTES)
    %let's use the default -> W
    MatTES=MaterialProperties('W');
end

if nargin < 7 | isempty(MatQET)
    %let's use the default -> Al
    MatQET=MaterialProperties('Al');
end

Q_TES= zeros(nTc,1);
Q_Rs = zeros(nTc,1);


Det=DetFunc(fridge,MatAbsb,MatTES,MatQET);
nRp=length(Det.elec.Rp);
Q_Rp =zeros(nTc,nRp);

for jj=1:nTc
    % let's change the Tc parameter in the W object:
    MatTES.Tc = Tc(jj);
    
    % Create the detector object:
    Det = DetFunc(fridge,MatAbsb,MatTES,MatQET);
    
    if Det.TES.Tc > Det.fridge.T_MC 
        % Now let's calculate the energy resolution:
        Det = Heating_Det(Det);
        
        Q_TES(jj) = Det.Heat.TES;
        Q_Rs(jj)  = Det.Heat.Rs;
        Q_Rp(jj,:)= Det.Heat.Rp;
    else
        Q_TES(jj) = NaN;
        Q_Rs(jj)  = NaN;
        Q_Rp(jj,:)= NaN;
    end    
end

%---------------- plot -------------------------------

cmap={'b','c','g','y'};

figure(1)
plot(Tc*1e3,Q_TES*1e12,'-k')
hold on
plot(Tc*1e3,Q_Rs*1e12,'-r')

leg_str={'TES: MC'
         ['Rs: ',Det.elec.Ts(3:end)]};
for jj=1:nRp
    leg_str=[leg_str;['Rp: ',Det.elec.Tp{jj}(3:end)]];
    plot(Tc*1e3,Q_Rp(:,jj)*1e12,['-',cmap{jj}])
end    
hold off
xlim([min(Tc),max(Tc)]*1e3)
grid on
set(gca,'yscale','log','xscale','log')
title([Det.name,': Fridge Heating vs TES Tc'])
ylabel('Fridge Heating [pW/Det]')
xlabel('TES Tc [mK]')
legend(leg_str,'location','northwest')

set(gca,'xtick',[10:10:100])



