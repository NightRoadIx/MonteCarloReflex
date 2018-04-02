%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [a,b] = FiberEllipseRadius (Tilt_Angle, Diameter)}
% Esta funci�n regresa los radios a & b de la elipse en la fibra.
% En el caso de que la fibra tenga una inclinaci�n cero, a = b
%
% Argumentos:
%   - Tilt_Angle = el �ngulo de inclinaci�n a partir de la normal de la
%   superficie del tejido. El valor debe estar dado en radianes.
%   - Diameter = El di�metro de la sonda en um
%
% Regresa:
%   - a: radio en la direcci�n X en um
%   - b: radio en la direcci�n Y en um
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