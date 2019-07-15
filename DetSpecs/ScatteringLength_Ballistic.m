function [ lscat ] = ScatteringLength_Ballistic(V,SA)
%
% In ballistic limited transport, the average scattering length will depend
% upon the absorber geometry
%
% From the super simple sabine equation (http://courses.physics.illinois.edu/phys406/Lecture_Notes/P406POM_Lecture_Notes/Derivation_of_the_Sabine_Equation.pdf)
%
% Inputs:
%           1) Absorber Volume
%           2) Absorber Surface Area
%
%
% 12/6/13: MCP
%--------------------------------------------------------------------------
lscat = 4*V/SA;
end

