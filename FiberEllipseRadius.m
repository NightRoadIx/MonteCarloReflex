%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [a,b] = FiberEllipseRadius (Tilt_Angle, Diameter)}
% Esta función regresa los radios a & b de la elipse en la fibra.
% En el caso de que la fibra tenga una inclinación cero, a = b
%
% Argumentos:
%   - Tilt_Angle = el ángulo de inclinación a partir de la normal de la
%   superficie del tejido. El valor debe estar dado en radianes.
%   - Diameter = El diámetro de la sonda en um
%
% Regresa:
%   - a: radio en la dirección X en um
%   - b: radio en la dirección Y en um
%
% Developed by: Roberto Reif / Boston University 
% Last Updated: September 3, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a,b] = FiberEllipseRadius(Tilt_Angle,Diameter)

% Calcula los radios de la elipse
d1 = cos(Tilt_Angle)*Diameter;
h = Diameter*sin(Tilt_Angle);
d2 = h/(tan(pi/2 - Tilt_Angle));
d3 = d1+d2;
    
a = d3/2;
b = Diameter/2;