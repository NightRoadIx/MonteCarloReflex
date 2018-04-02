%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BHMIE
% Notes: Basado en Bohren y Huffman Apéndice A
% Desarrollado por: Roberto Reif
% Fecha: 05/18/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula los elementos de amplitud de la matriz de esparcimiento y
% eficiencias para la extinción, esparcimiento total y retroesparcimiento
% para un parámetro de tamaño dado e índice de refracción relativo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S1,S2] = BHMIE(X,REFREL,NANG)

% Decidir cual será el valor de NMX, explicado en el apéndice A
DX = X;
Y = X*REFREL;
XSTOP = (X+4^(1/3)+2);
NSTOP = XSTOP;
YMOD = abs(Y);

NMX = ceil(max([XSTOP YMOD]) + 15);
    
DANG = pi/(2*(NANG - 1));

% Obtiene todos los ángulo y su coseno
THETA = [];

for J=1:NANG 
    THETA(J) = (J-1)*DANG;
    AMU(J) = cos(THETA(J));
end


% Calcula la derivada logaritmica D(J) por recurrencia "hacia abajo"
% comenzando con el valor inicial de 0 en J = NMX
D = [];
D(NMX) = 0.0001 + 0.0001*i;
NN = NMX - 1;
for N = 1:NN
    RN = NMX - N + 1;
    D(NMX - N) = (RN/Y)-(1/(D(NMX - N + 1) + (RN / Y)));
end

PI0 = [];
PI1 = [];
S1 = [];
S2 = [];
for J = 1:NANG
    PI0(J) = 0;
    PI1(J) = 1;
end

NN = 2*NANG - 1;
for J = 1:NN
    S1(J) = 0 + 0i;
    S2(J) = 0 + 0i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCIONES DE RICATTI-BESSEL CON ARGUMENTO REAL X CALCULADO POR
% RECURRENCIA "HACIA ARRIBA"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI0 = cos(DX);
PSI1 = sin(DX);
CHI0 = -sin(X);
CHI1 = cos(X);
APSI0 = PSI0;
APSI1 = PSI1;
XI0 = APSI0 - CHI0*i;
XI1 = APSI1 - CHI1*i;
QSCA = 0;
N = 1;

AN = [];
BN = [];

for N=1:NSTOP
    DN = N;
    RN = N;
    FN = (2*RN + 1)/(RN * (RN+1));
    PSI = (2*DN - 1) * PSI1/DX - PSI0;
    APSI = PSI;
        CHI = (2 * RN - 1) * CHI1/X - CHI0;
    XI = APSI - CHI * i;
    AN(N) = ((D(N)/REFREL) + RN/X)*APSI - APSI1;
    AN(N) = AN(N)/(((D(N)/REFREL) + RN/X)*XI - XI1);
    BN(N) = ((REFREL * D(N)) + RN /X)*APSI - APSI1;
    BN(N) = BN(N)/(((REFREL*D(N))+RN/X)*XI - XI1);
    QSCA = QSCA + (2*RN+1)*((abs(AN(N)))^2 + (abs(BN(N)))^2);

    for J=1:NANG
        JJ = 2*NANG - J;
        PI(J) = PI1(J);
        TAU(J) = RN*AMU(J) * PI(J) - (RN+1)*PI0(J);
        P = (-1)^(N-1);
        S1(J) = S1(J) + FN*(AN(N)*PI(J)+BN(N)*TAU(J));
        T = (-1)^N;
        S2(J) = S2(J) + FN * (AN(N)*TAU(J) + BN(N)*PI(J));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if J~=JJ
            S1(JJ) = S1(JJ) + FN * (AN(N)*PI(J)*P+BN(N)*TAU(J)*T);
            S2(JJ) = S2(JJ) + FN*(AN(N)*TAU(J)*T+BN(N)*PI(J)*P);
        end
    end

    PSI0 = PSI1;
    PSI1 = PSI;
    APSI1 = PSI1;
    CHI0 = CHI1;
    CHI1 = CHI;
    XI1 = APSI1 - CHI1*i;
    RN = N + 1;

    for J =1:NANG
        PI1(J) = ((2*RN-1)/(RN-1))*AMU(J)*PI(J);
        PI1(J) = PI1(J) - RN*PI0(J)/(RN - 1);
        PI0(J) = PI(J);
    end

end

QSCA = (2/(X*X))*QSCA;
QEXT = (4/(X*X))*real(S1(1));
QBACK = (4/(X*X))*abs(S1(2*NANG-1))*abs(S1(2*NANG-1));
