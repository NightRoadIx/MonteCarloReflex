%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ux,uy,uz,tita] = MIEaprox(ux,uy,uz,MieAng,MiePhase);
% This function returns the ux, uy and uz direction of the photon once it
% hit a scatterer. The directional vector is referenced with XYZ of the medium. 
%
% Arguments:
%    - ux: direction in x
%    - uy: direction in y
%    - uz: direction in z
%    - g: the anisotropy factor
%
% Returns:
%    - ux: direction in x
%    - uy: direction in y
%    - uz: direction in z
%    - tita: scattering angle
%
% Developed by: Roberto Reif / Boston University 
% Last Updated: September 07, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ux,uy,uz,tita] = MIEaprox (ux,uy,uz,MieAng,MiePhase)

% adds the Mie Phase, the last number should be 1
Prob = cumsum(MiePhase);
x = rand;

% gets the value where Prob is equal or greater than the random number
LocVal = find(Prob>=x);

% gets the first value where Prob >= x and that is the angle that is
% selected
tita = MieAng(LocVal(1))*pi/180;

Azim_Angle = 2*pi*rand;

temp_ux = ux;
temp_uy = uy;
temp_uz = uz;

% this is taken from page 15 of "MonteCarlo Modeling of Light Transport in MultiLayered Tissue in Standard C" Wang-Jacques, 
% I need to check this

if uz ~= 1
    ux = sin(tita)/sqrt(1-temp_uz^2)*(temp_ux*temp_uz*cos(Azim_Angle) - temp_uy*sin(Azim_Angle)) + temp_ux*cos(tita);
    uy = sin(tita)/sqrt(1-temp_uz^2)*(temp_uy*temp_uz*cos(Azim_Angle) + temp_ux*sin(Azim_Angle)) + temp_uy*cos(tita);
    uz = - sin(tita)*cos(Azim_Angle)*sqrt(1-temp_uz^2) + temp_uz*cos(tita);
else
    ux = sin(tita)*cos(Azim_Angle);
    uy = sin(tita)*sin(Azim_Angle);
    uz = sign(temp_uz)*cos(tita);
end

