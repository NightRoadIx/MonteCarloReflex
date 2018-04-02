clear all
close all
clc

tic
% Tomar el tiempo como la "semilla" para iniciar el generador de números
% aleatorios
rand('state',sum(100*clock));

% CONSTANTES
YES = 1;
NO = 0;
% Definir al constante de Planck
hplanck = 6.626068*10^-34;
% Velocidad de la luz en el vacío
cluz = 299792.458; % m/s

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% PROPIEADES DEL TEJIDO
% Longitud de onda de análisis
lambda = 633;
% Coeficiente de esparcimiento
us = 143.7;    % (1/cm)
% Coeficiente de absorción
ua = 5.9; % (1/cm)
% Apertura numérica de la fibra
NA = 0.22;
% Indice de refracción de la fibra
n_fiber = 1.42;
% Indice de refracción del medio
n_medium = 1.37; % water = 1.37
% Velocidad de la luz en el medio
cmedium = cluz/n_medium;
% Paso fijo
Fix_Step = NO;
% Utilizar ecuación Henyey-Greenstein
HG = YES;

if HG == YES % En caso de que sea x aproximación de HG
    g = 0.95;       %Anisotropía
    MieAng = 0;     %Parámetro ángulo 
    MiePhase = 0;   %Parámetro fase 
else % En el caso de Mie
    n_sphere = 1.33; % Esfera. Esto puede ser un número complejo
    sphere_radius = 0.26; %Radio de la esfera en um
    wavel = lambda/1000; % Longitud de onda en um
    % Calcular el parámetro de tamaño "X" y el indice de refracción relativo
    % para una esfera de índice de refracción conocido, índice de refracción
    % del medio conocido, radio y longitud de onda en el espacio...     
    [MieAng,MiePhase,g] = CALLBH(n_medium,real(n_sphere),imag(n_sphere),sphere_radius,wavel);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Profunidad de penetración óptica
ydelta = 1/sqrt( 3*ua*(ua+us*(1-g)) );

Zfibers = 0; % El valor de esto debe de ser >= 0 (Eje z)

Source_Tilt_Angle = 0*pi/180;   % ángulo de inclinación de la fuente en radianes; negativo hacia el detector
Detector_Tilt_Angle = 0*pi/180; % ángulo de inclinación del detector en radianes; positivo a la fuente

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% ESTRUCUTURA DE LA FIBRA DE PRUEBA
Num_Sources = 6;            % Número de fuentes
Source_Diameter = 400;      % Diámetro de la fuente en um
Detector_Diameter = 400;    % Diámetro del detector en um
Fiber_Separation = 420;     % Separación de las fibras en um
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Esta es la máxima separación que existe
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
    Fiber(n).SizeCoord = 1;     % Tamaño de coordenadas
    Fiber(n).PhotonCounter = 0; % Contador de fotones
    Fiber(n).SizeScatAngle = 1; % Ángulo de esparcimiento 
end

% %Comenzar la simulación
Xfin = [];
Yfin = [];
Zfin = [];

%Obtener el radio de la fuente y del detector
source_radius = Source_Diameter/2;
detect_radius = Detector_Diameter/2;

% Obtener la función de lo(s) detector(es)
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
        
    % Obtener la función de la(s) fuente(s)
    th = 0:pi/50:2*pi;
    xunit = source_radius * cos(th) + x(k);
    yunit = source_radius * sin(th) + y(k);
    
    % Graficar la(s) fuente(s)
    plot3(xunit,yunit,yunit-yunit,'color', A(k,:))
    
end
% Crear un arreglo con los centros de las fuentes
cSource = [x', y'];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Número de fotones a recolectar por fuente
NumPhotonsCollected = 500;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % COMENZAR LA SIMULACIÓN % % % % % % % 
% Lanzar los fotones
for counter = 1: NumPhotonsCollected
    for qu = 1:length(cSource)
       [XSave,YSave,ZSave,PhotonsLaunched, ScatAngle, fiber] = PropagatePhotonb(us,ua,Fix_Step,HG,g,aSource,bSource,cSource(qu,:),aDetector,bDetector,Zfibers, Fiber_Separation, NA, Source_Tilt_Angle, Detector_Tilt_Angle, n_fiber, n_medium, MaxDistance, MieAng, MiePhase);
       
       % Sellar los valores del movimiento de cada fotón registrado
       % Xfin=[XSave];
       % Yfin=[YSave];
       % Zfin=[ZSave];
       
       % Incrementar el contador de fotón
       Fiber(fiber).PhotonCounter = Fiber(fiber).PhotonCounter + 1;
       % Posicion del fotón
       Fiber(fiber).PhotonPosition(Fiber(fiber).PhotonCounter) = Fiber(fiber).SizeCoord;
       
       % Grabar la posición XYZ
       for n = Fiber(fiber).SizeCoord:Fiber(fiber).SizeCoord+size(XSave,2)-1
           Fiber(fiber).XYZ(n,:) = [XSave(n-Fiber(fiber).SizeCoord+1) YSave(n-Fiber(fiber).SizeCoord+1) ZSave(n-Fiber(fiber).SizeCoord+1)];
       end
       
       % Grabar los ángulos
       for n = Fiber(fiber).SizeScatAngle:Fiber(fiber).SizeScatAngle+size(ScatAngle,2)-1
           Fiber(fiber).ScatAngle(n,:) = ScatAngle(n-Fiber(fiber).SizeScatAngle+1);
       end
       
       Fiber(fiber).SizeCoord = Fiber(fiber).SizeCoord+size(XSave,2);
       Fiber(fiber).SizeScatAngle = Fiber(fiber).SizeScatAngle+size(ScatAngle,2);
       %Incrementar el contador de fotones lanzados
       NumPhotonsLaunched = NumPhotonsLaunched+PhotonsLaunched;
       disp(strcat('Fuente: ', int2str(qu)));
       disp(strcat('# de fotones lanzados: ',int2str(NumPhotonsLaunched)))
       
       % Si ha habido un número de fotones lanzados entonces se termina el
       % código (Nmax = 30 000 000)
       if NumPhotonsLaunched > 3000000000
          break
      end
      disp(strcat('Fotón Capturado ',int2str(counter)))
      
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

%Crear los títulos de las gráficas
xlabel('x [\mum]');
ylabel('y [\mum]');
zlabel('z [\mum]');
cad1='Simulación Monte Carlo reflectancia difusa';
cad2=strcat('D_f_u_e_n_t_e = ',num2str(Source_Diameter),'mm, A_d_e_t = 7 mm^2, dist_D_-_E =  ',num2str(Fiber_Separation));
cad3=strcat('\mu_a =  ',num2str(ua),' mm^-^1','\mu_s =  ',num2str(us),' mm^-^1, g =  ',num2str(g),', n =  ',num2str(n_medium));
cad4=strcat('N_f_o_t_o_n_e_s =  ', num2str(Fiber.PhotonCounter));
title(strcat(cad1,'\newline',cad2,'\newline',cad3,'\newline',cad4),'FontWeight','bold')
grid

%Terminar de contar el tiempo de simulación
T1=toc;
disp(' '), disp(cad1), disp(cad2), disp(cad3), disp(cad4)
disp(strcat('Tiempo simulación: ',' ',num2str(T1),' segundos'))

%Salvar los datos simulados
save( strcat('Proba',num2str(lambda)) )