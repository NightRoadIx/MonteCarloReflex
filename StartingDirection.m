%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ux,uy,uz] = StartingDirection (NA, Tilt_Angle, n_fiber, n_medium)
% This function returns the ux, uy and uz direction of the photon exiting
% the fiber probe
%
% Arguments:    
%    - NA = Numerical Aperture or sin(teta) of the fiber referenced to air.
%    - Tilt_Angle = Angle of tilt of the fiber in radians.
%               positive -> away from detector, negative -> towards detector
%    - n_fiber = index of refraction of the fiber
%    - n_medium = index of refraction of the medium
%
% Returns:
%   - ux: X direction
%   - uy: Y direction
%   - uz: Z direction
% If there is an invalid combination of the input arguments, then the
% return will be ux, uy, uz = 10!
%
% Developed by: Roberto Reif / Boston University 
% Last Updated: June 21, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ux,uy,uz] = StartingDirection (NA, Tilt_Angle, n_fiber, n_medium)

% gets the critical angle inside the fiber
tita_fiber = asin(NA/n_fiber);

% gets the angles in the y direction
titay1 = asin(NA/n_medium);
titay2 = -titay1;

% gets the angles in the x direction
titax1 = asin(n_fiber/n_medium*sin(-Tilt_Angle+tita_fiber));
titax2 = asin(n_fiber/n_medium*sin(-Tilt_Angle-tita_fiber));

if imag(titax1) ~= 0|| imag (titax2) ~= 0
    input('Error Starting Direction')
    return
end

titaNA = (titax1 - titax2)/2;
titaHalf = (titax1 + titax2)/2;

a = sin(titaNA);
b = sin(titay1);

ux = rand*a*sign(randn);
uy = rand*b*sign(randn)*sqrt(1-(ux/a)^2);
uz = sqrt(1-ux^2-uy^2);

xztita = atan(ux/(uz+eps));
h = uz/(cos(xztita)+eps);
ux = sin(xztita+titaHalf)*h;
uz = cos(xztita+titaHalf)*h;
