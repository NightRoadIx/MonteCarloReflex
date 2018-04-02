clear all
close all
clc

tic
% Tomar el tiempo como la "semilla" para iniciar el generador de n�meros
% aleatorios
rand('state',sum(100*clock));

% CONSTANTES
YES = 1;
NO = 0;
% Definir al constante de Planck
hplanck = 6.626068*10^-34;
% Velocidad de la luz en el vac�o
cluz = 299792.458; % m/s

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% PROPIEADES DEL TEJIDO
% Longitud de onda de an�lisis
lambda = 633;
% Coeficiente de esparcimiento
us = 143.7;    % (1/cm)
% Coeficiente de absorci�n
ua = 5.9; % (1/cm)
% Apertura num�rica de la fibra
NA = 0.22;
% Indice de refracci�n de la fibra
n_fiber = 1.42;
% Indice de refracci�n del medio
n_medium = 1.37; % water = 1.37
% Velocidad de la luz en el medio
cmedium = cluz/n_medium;
% Paso fijo
Fix_Step = NO;
% Utilizar ecuaci�n Henyey-Greenstein
HG = YES;

if HG == YES % En caso de que sea x aproximaci�n de HG
    g = 0.95;       %Anisotrop�a
    MieAng = 0;     %Par�metro �ngulo 
    MiePhase = 0;   %Par�metro fase 
else % En el caso de Mie
    n_sphere = 1.33; % Esfera. Esto puede ser un n�mero complejo
    sphere_radius = 0.26; %Radio de la esfera en um
    wavel = lambda/1000; % Longitud de onda en um
    % Calcular el par�metro de tama�o "X" y el indice de refracci�n relativo
    % para una esfera de �ndice de refracci�n conocido, �ndice de refracci�n
    % del medio conocido, radio y longitud de onda en el espacio...     
    [MieAng,MiePhase,g] = CALLBH(n_medium,real(n_sphere),imag(n_sphere),sphere_radius,wavel);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Profunidad de penetraci�n �ptica
ydelta = 1/sqrt( 3*ua*(ua+us*(1-g)) );

Zfibers = 0; % El valor de esto debe de ser >= 0 (Eje z)

Source_Tilt_Angle = 0*pi/180;   % �ngulo de inclinaci�n de la fuente en radianes; negativo hacia el detector
Detector_Tilt_Angle = 0*pi/180; % �ngulo de inclinaci�n del detector en radianes; positivo a la fuente

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% ESTRUCUTURA DE LA FIBRA DE PRUEBA
Num_Sources = 6;            % N�mero de fuentes
Source_Diameter = 400;      % Di�metro de la fuente en um
Detector_Diameter = 400;    % Di�metro del detector en um
Fiber_Separation = 420;     % Separaci�n de las fibras en um
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Esta es la m�xima separaci�n que existe
MaxDistance = max(4*max(Fiber_Separation),Source_Diameter);

% Obtener el radio de la elipse en la punta de la fibra
% Fuente(s)
[aSource,bSource] = FiberEllipseRadius(Source_Tilt_Angle,Source_Diameter);
% Detector
[aDetector,bDetector] = FiberEllipseRadius(Detector_Tilt_Angle,Detector_Diameter);

NumPhotonsLaunched = 0;     % iniciar el contador de # de fotones lanzados

%
% Llenar los vectores de la estructura que representa a la fibra
% (esto es para cuando hay varios detectores)
for n = 1:size(Fiber_Separation,2)
    Fiber(n).SizeCoord = 1;     % Tama�o de coordenadas
    Fiber(n).PhotonCounter = 0; % Contador de fotones
    Fiber(n).SizeScatAngle = 1; % �ngulo de esparcimiento 
end

% %Comenzar la simulaci�n
Xfin = [];
Yfin = [];
Zfin = [];

%Obtener el radio de la fuente y del detector
source_radius = Source_Diameter/2;
detect_radius = Detector_Diameter/2;

% Obtener la funci�n de lo(s) detector(es)
t2 = ( (Fiber_Separation/2) - detect_radius ):(detect_radius + (Fiber_Separation/2));
circ2 = sqrt(detect_radius^2 - (t2 - (Fiber_Separation/2) ).^2);
% Graficar el detector
plot3(t2,circ2,circ2-circ2,'k',t2,-circ2,circ2-circ2,'k')
hold on; grid

% Centro del detector
x1 = Fiber_Separation/2;
y1 = 0;
x = zeros(1,Num_Sources);
y = zeros(1,Num_Sources);
A = hsv(Num_Sources);
for k=1:Num_Sources
    % Obtener los centros de las fuentes
    x(k) = fix(x1 + Fiber_Separation*sin( (2/Num_Sources)*pi*(k-1) ));
    y(k) = fix(y1 + Fiber_Separation*cos( (2/Num_Sources)*pi*(k-1) ));
        
    % Obtener la funci�n de la(s) fuente(s)
    th = 0:pi/50:2*pi;
    xunit = source_radius * cos(th) + x(k);
    yunit = source_radius * sin(th) + y(k);
    
    % Graficar la(s) fuente(s)
    plot3(xunit,yunit,yunit-yunit,'color', A(k,:))
    
end
% Crear un arreglo con los centros de las fuentes
cSource = [x', y'];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% N�mero de fotones a recolectar por fuente
NumPhotonsCollected = 500;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % COMENZAR LA SIMULACI�N % % % % % % % 
% Lanzar los fotones
for counter = 1: NumPhotonsCollected
    for qu = 1:length(cSource)
       [XSave,YSave,ZSave,PhotonsLaunched, ScatAngle, fiber] = PropagatePhotonb(us,ua,Fix_Step,HG,g,aSource,bSource,cSource(qu,:),aDetector,bDetector,Zfibers, Fiber_Separation, NA, Source_Tilt_Angle, Detector_Tilt_Angle, n_fiber, n_medium, MaxDistance, MieAng, MiePhase);
       
       % Sellar los valores del movimiento de cada fot�n registrado
       % Xfin=[XSave];
       % Yfin=[YSave];
       % Zfin=[ZSave];
       
       % Incrementar el contador de fot�n
       Fiber(fiber).PhotonCounter = Fiber(fiber).PhotonCounter + 1;
       % Posicion del fot�n
       Fiber(fiber).PhotonPosition(Fiber(fiber).PhotonCounter) = Fiber(fiber).SizeCoord;
       
       % Grabar la posici�n XYZ
       for n = Fiber(fiber).SizeCoord:Fiber(fiber).SizeCoord+size(XSave,2)-1
           Fiber(fiber).XYZ(n,:) = [XSave(n-Fiber(fiber).SizeCoord+1) YSave(n-Fiber(fiber).SizeCoord+1) ZSave(n-Fiber(fiber).SizeCoord+1)];
       end
       
       % Grabar los �ngulos
       for n = Fiber(fiber).SizeScatAngle:Fiber(fiber).SizeScatAngle+size(ScatAngle,2)-1
           Fiber(fiber).ScatAngle(n,:) = ScatAngle(n-Fiber(fiber).SizeScatAngle+1);
       end
       
       Fiber(fiber).SizeCoord = Fiber(fiber).SizeCoord+size(XSave,2);
       Fiber(fiber).SizeScatAngle = Fiber(fiber).SizeScatAngle+size(ScatAngle,2);
       %Incrementar el contador de fotones lanzados
       NumPhotonsLaunched = NumPhotonsLaunched+PhotonsLaunched;
       disp(strcat('Fuente: ', int2str(qu)));
       disp(strcat('# de fotones lanzados: ',int2str(NumPhotonsLaunched)))
       
       % Si ha habido un n�mero de fotones lanzados entonces se termina el
       % c�digo (Nmax = 30 000 000)
       if NumPhotonsLaunched > 3000000000
          break
      end
      disp(strcat('Fot�n Capturado ',int2str(counter)))
      
      plot3(XSave,YSave,ZSave,'.', 'color', A(qu,:));
      drawnow;
%       plot3(200-XSave,YSave,ZSave,'r.-')
%       plot3(XSave,YSave,ZSave,XSave,YSave,ZSave,'o')
      hold on
    end
end

for f = 1:size(Fiber_Separation,2)
    if Fiber(f).PhotonCounter~=0
        Fiber(f).PhotonPosition(Fiber(f).PhotonCounter+1) = Fiber(f).SizeCoord;
    end
end

%Crear los t�tulos de las gr�ficas
xlabel('x [\mum]');
ylabel('y [\mum]');
zlabel('z [\mum]');
cad1='Simulaci�n Monte Carlo reflectancia difusa';
cad2=strcat('D_f_u_e_n_t_e = ',num2str(Source_Diameter),'mm, A_d_e_t = 7 mm^2, dist_D_-_E =  ',num2str(Fiber_Separation));
cad3=strcat('\mu_a =  ',num2str(ua),' mm^-^1','\mu_s =  ',num2str(us),' mm^-^1, g =  ',num2str(g),', n =  ',num2str(n_medium));
cad4=strcat('N_f_o_t_o_n_e_s =  ', num2str(Fiber.PhotonCounter));
title(strcat(cad1,'\newline',cad2,'\newline',cad3,'\newline',cad4),'FontWeight','bold')
grid

%Terminar de contar el tiempo de simulaci�n
T1=toc;
disp(' '), disp(cad1), disp(cad2), disp(cad3), disp(cad4)
disp(strcat('Tiempo simulaci�n: ',' ',num2str(T1),' segundos'))

%Salvar los datos simulados
save( strcat('Proba',num2str(lambda)) )