%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBH
% Notas: Basado en Bohren y Huffman Apéndice A
% Desarrollado por: Roberto Reif
% Fecha: 05/18/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula el parámetro de tamaño "X" y el índice de refracción relativo
% (REFREL) para un indice de refracción de esfera, índice de refracción del
% medio, radio y longitud de onda en espacio libre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ctes
% REFMED = n del medio
% REFRE = Parte real de n de la esfera
% REFIM = Parte Imaginaria de n de la esfera
% RAD = Radio de la esfera (um)
% WAVEL = Longitud de onda (um)
% NANG = Número de ángulos entre 0-90 grados
% DANG = 90(rad)/(NANG-1) = Radianes por división

function [ANG,NOTPOL,G] = CALLBH(REFMED,REFRE,REFIM,RAD,WAVEL)

NANG = 500;
DANG = pi/(2*(NANG-1));

% REFREL = índice de refracción relativo
% X = Parámetro del tamaño
REFREL = (REFRE+REFIM*i)/REFMED;
X = 2*pi*RAD*REFMED/WAVEL;

[S1,S2] = BHMIE(X,REFREL,NANG);

S11NORM = 0.5*(abs(S2).^2 +abs(S1).^2);
S11NOR = sum(S11NORM);
NAN = 2*NANG - 1;
for J = 1:NAN
    S11(J) = 0.5*abs(S2(J))^2;
    S11(J) = S11(J) + 0.5 *abs(S1(J))^2;
    S12(J) = 0.5*abs(S2(J))^2;
    S12(J) = S12(J) - 0.5*abs(S1(J))^2;
    POL(J) = -S12(J)/S11(J);
    S33(J) = real(S2(J)*conj(S1(J)));
    S33(J) = S33(J)/S11(J);
    S34(J) = imag(S2(J)*conj(S1(J)));
    S34(J) = S34(J)/S11(J);
    S11(J) = S11(J)/S11NOR;
    S12(J) = S12(J)/S11NOR;%%%%
    S33(J) = S33(J)/S11NOR;%%%%
    S34(J) = S34(J)/S11NOR;%%%%
    ANG(J) = DANG*(J-1)*180/pi;
end

IPAR = S11 + S12;
IPER = S11 - S12;
NOTPOL = S11;
P = -S12./S11;

RADANG = ANG.*pi/180;

G = 0;
GNORM = 0;
for I = 1:NAN
   G = G + cos(RADANG(I))*NOTPOL(I);
   GNORM = GNORM + NOTPOL(I);
end
G = G./GNORM;