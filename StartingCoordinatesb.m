%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,Y] = StartingCoordinates (a,b,Fiber_Separation)
% Esta función regresa las coordenadas X & Y del fotón saliendo de la
% fuente de la fibra
%
% Argumentos:
%    - a, b = radios de la elipse
%    - c = punto central de la fuente
%    - Fiber_Separation =  la separación entre la fuente y el detector (um)
%
% Returns:
%   - X: la coordenada en X
%   - Y: la coordenada en Y
%
% Developed by: Roberto Reif / Boston University 
% Last Updated: July 19, 2004
%               July  9, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y] = StartingCoordinatesb(a,b,c,Fiber_Separation)

% % inicia X, Y a un valor inválido
% X = a^2-Fiber_Separation^2;
% Y = b^2;

X = c(1) + a*(2*rand-1);
Y = sqrt((1 - (X - c(1))^2/a^2)*b^2)*(2*rand-1) + c(2);

% switch c
%     case 1
%         % Fibra fuente localizada en la parte negativa del eje x
%         X = a*(2*rand-1) - Fiber_Separation/2;
%         Y = sqrt((1 - (X + Fiber_Separation/2)^2/a^2)*b^2)*(2*rand-1);        
%     case 2
%         % Fibra fuente localizada en la parte positiva del eje y
%         X = a*(2*rand-1) + Fiber_Separation/2;
%         Y = sqrt((1 - (X - Fiber_Separation/2)^2/a^2)*b^2)*(2*rand-1) + Fiber_Separation;        
%     case 3
%         % Fibra fuente localizada en la parte positiva del eje x
%         X = a*(2*rand-1) + (3/2)*Fiber_Separation;
%         Y = sqrt((1 - (X - (3/2)*Fiber_Separation)^2/a^2)*b^2)*(2*rand-1);        
%     case 4
%         % Fibra fuente localizada en la parte negativa del eje y
%         X = a*(2*rand-1) + Fiber_Separation/2;
%         Y = sqrt((1 - (X - Fiber_Separation/2)^2/a^2)*b^2)*(2*rand-1) - Fiber_Separation;        
% end
% % Fibra fuente localizada en la parte negativa del eje x
% X = a*(2*rand-1) - Fiber_Separation/2;
% Y = sqrt((1 - (X + Fiber_Separation/2)^2/a^2)*b^2)*(2*rand-1);

% % Fibra fuente localizada en la parte positiva del eje x
% X = a*(2*rand-1) + (3/2)*Fiber_Separation;
% Y = sqrt((1 - (X + Fiber_Separation/2)^2/a^2)*b^2)*(2*rand-1);
% 
% % Fibra fuente localizada en la parte negativa del eje y
% X = a*(2*rand-1) - Fiber_Separation/2;
% Y = -Fiber_Separation + sqrt((1 - (X + Fiber_Separation/2)^2/a^2)*b^2)*(2*rand-1);
% 
% % Fibra fuente localizada en la parte positiva del eje y
% X = a*(2*rand-1) - Fiber_Separation/2;
% Y = Fiber_Separation + sqrt((1 - (X + Fiber_Separation/2)^2/a^2)*b^2)*(2*rand-1);
