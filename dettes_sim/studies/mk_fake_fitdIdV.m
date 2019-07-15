function [fitdIdV]= mk_fake_fitdIdV(Det)
%this function makes a fake dIdV fit to study properties of the degeneracy

fitdIdV=[];

fitdIdV.dIdV0 = Det.Response.dIdV0;
fitdIdV.w_I   = Det.TES.w_I;
fitdIdV.wp_p  = Det.TES.wp_p;
fitdIdV.wp_m  = Det.TES.wp_m;
fitdIdV.Rl    = Det.elec.Rl;
fitdIdV.Ro    = Det.TES.Ro;
end