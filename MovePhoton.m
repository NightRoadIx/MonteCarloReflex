%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,Y,Z] = MovePhoton(X, Y, Z, ux, uy, uz, Fix_Step, us, ua);
% Esta función regresa las coordenadas X, Y & Z del fotón una vez que este
% ha dado un paso para propagarse
%
% Argumentos:
%    - X: posición original X (um)
%    - Y: posición original Y (um)
%    - Z: posición original Z (um)
%    - ux: dirección en x
%    - uy: dirección en y
%    - uz: dirección en z
%    - Fix_Step: 
%               1 si el fotón se mueve 1/ut 
%               0 si el fotón se mueve un promedio de 1/ut
%    - us: Coeficiente de esparcimiento (1/cm)
%    - ua: Coeficiente de absorción (1/cm)
%
% Returns:
%   - X: la coordenada X (in um)
%   - Y: la coordenada Y (in um)
%   - Z: la coordenada Z (in um)
%
% Desarrollado por: Roberto Reif / Boston University 
% Última actualización: Junio 29, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X1,Y1,Z1] = MovePhoton(X, Y, Z, ux, uy, uz, Fix_Step, us, ua, Zfibers)

    [step] = Stepsize (us,ua, Fix_Step);
    
    % Calcular el incremento de posición
    dx = ux*step; 
    dy = uy*step;
    dz = uz*step;

    % Calcular la nueva posición
    X1 = X+dx; 
    Y1 = Y+dy;
    Z1 = Z+dz;
    
    if Z1<Zfibers
        Proportion = Z/(Z-Z1);
        X1=X-(X-X1)*Proportion;
        Y1=Y-(Y-Y1)*Proportion;
        Z1=Zfibers;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DS] = Stepsize (us,ua)
% Esta función regresa el tamaño del paso del fotón hasta que este
% interactúa con un esparcidor o un absorbente
%
% Argumentos:    
%    - us = Coeficiente de esparcimiento (cm-1)
%    - ua = Coeficiente de absorción (cm-1)
%        
% Regresa:
%   - DS: Tamaño del paso en um
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DS] = Stepsize (us,ua,Fixed)

YES = 1;

if Fixed == YES
    DS = 1/(us+ua)*10000; % Regresa en um
else
    DS = -log(rand+eps)/(us+ua)*10000; % Esto es en ln
end