%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [value] = DetectDirection (NA, Tilt_Angle, n_fiber, n_medium, uxout, uyout, uzout)
% Esta funci�n regresa un SI si el fot�n est� dentro de NA de la fibra, en
% otro caso regresa NO
%
% Argumentos:    
%    - NA = Apertura num�rica o sen(theta) de la fibra referenciado al aire
%    - Tilt_Angle = �ngulo de inclinaci�n de la fibra en radianes.
%               positivo -> alej�ndose del detector 
%               negativo -> hacia el detector
%    - n_fiber = �ndice de refracci�n de la fibra
%    - n_medium = �ndice de refracci�n del medio
%    - uxout, uyout, uzout = Direcci�n del fot�n en x, y & z
%
% Regresa:
%   - value = 1 si el fot�n esta dentro de la fibra, 0 en otro caso
%
% Desarrolado por: Roberto Reif / Boston University 
% �ltima actualizaci�n: Septiembre 07, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value] = DetectDirection (NA, Tilt_Angle, n_fiber, n_medium,uxout,uyout,uzout)

YES=1;
NO=0;

% Obtiene el �ngulo cr�tico dentro de la fibra
tita_fiber = asin(NA/n_fiber);

% Obtiene los �ngulos en la direcci�n y
titay1 = asin(NA/n_medium);
titay2 = -titay1;

% Obtiene los �ngulos en la direcci�n x
titax1 = asin(n_fiber/n_medium*sin(-Tilt_Angle+tita_fiber));
titax2 = asin(n_fiber/n_medium*sin(-Tilt_Angle-tita_fiber));

titaNA = (titax1 - titax2)/2;
titaHalf = (titax1 + titax2)/2;

a = sin(titaNA);
b = sin(titay1);

% Obtiene los �ngulos en la direcci�n YZ
yangle = atan(uyout/(uzout+eps));
uy = uyout;

% Revisa si el fot�n est� dentro del �ngulo YZ
if yangle<=titay1 && yangle>=titay2
    
    % Obtiene el �ngulo en la direcci�n XZ
    ux1 = a*sqrt(1-(uy/b)^2);
    uz1 = sqrt(1-ux1^2-uy^2);
    
    xztita = atan(ux1/(uz1+eps));
    h = uz1/(cos(xztita)+eps);
    ux1 = sin(xztita+titaHalf)*h;
    ux2 = sin(-xztita+titaHalf)*h;
    uz1 = cos(xztita+titaHalf)*h;
    uz2 = cos(-xztita+titaHalf)*h;
    
    titax1 = atan(ux1/uz1);
    titax2 = atan(ux2/uz2);
    xangle = atan(uxout/(uzout+eps));
    
    % Verifica si el �ngulo XZ es correcto
    if (xangle<=titax1 && xangle>=titax2) || (xangle>=titax1 && xangle<=titax2)
        value = YES;
        return
    end    
end

value = NO;
