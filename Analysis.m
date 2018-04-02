close all
clear all
clc

% Recargar los resultados
uiload

% Definir el número de Voxels buscados
Number_Voxels = 100;
%Definir al constante de Planck
hplanck = 6.626068*10^-34;

%Profunidad de penetración óptica
ydelta = 1/sqrt( 3*ua*(ua+us*(1-g)) );

%Obtener el radio de la fuente y del detector
source_radius=Source_Diameter/2;
detect_radius=Detector_Diameter/2;

%t1 = -Source_Diameter:-0;
t1 = -(source_radius + (Fiber_Separation/2)):-( (Fiber_Separation/2) - source_radius );
circ1 = sqrt(source_radius^2 - (t1 + (Fiber_Separation/2) ).^2);
%t2 = 0:Source_Diameter;
t2 = ( (Fiber_Separation/2) - detect_radius ):(detect_radius + (Fiber_Separation/2));
circ2 = sqrt(detect_radius^2 - (t2 - (Fiber_Separation/2) ).^2);
%t3=Source_Diameter:2*Source_Diameter;

%Graficar la circunferencia que representa al detector
figure
plot3(t2,circ2,circ2-circ2, 'k',t2,-circ2,circ2-circ2,'k'), grid
hold on

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
%     A(k,:) = 0.5 + 0.5.*[rand rand rand];
    plot3(xunit,yunit,yunit-yunit,'color', A(k,:))
    
end

% Sobre la estrcutura con los datos de la simulación
for fibercounter = 1:size(Fiber,2)
    % MaxDistance = max(4*max(Fiber_Separation),2000);
    
    %Si el número de fotones en el detector es mayor a 0
    if Fiber(fibercounter).PhotonCounter>0

        % Se obtiene la distancia máxima de separación de fuente-detector
        MaxDistance = max(4*(Fiber_Separation(fibercounter)),2000);
        %Tamaño del Voxel
        Voxel_Size(fibercounter) = 2*ceil (MaxDistance / Number_Voxels);

        %Limpiar Fiber.Matrix
        clear Fiber(fibercounter).Matrix
        %Rellenar esta Matrix con ceros
        Fiber(fibercounter).Matrix(Number_Voxels,Number_Voxels,Number_Voxels) = zeros;

        %Ir recorriendo todos los fotones simulados que llegaron al
        %detector
        disp('Fotones simulados'), k = 1;
        for i = 1:Fiber(fibercounter).PhotonCounter
            disp(strcat('Foton  ',num2str(i)))
            clear XS1, clear YS1, clear ZS1
            % Obtener las coordenadas XYZ del camino del fotón
            for f = Fiber(fibercounter).PhotonPosition(i):Fiber(fibercounter).PhotonPosition(i+1)-1
%                 disp(strcat('   ->',num2str(f)));
                XS1(f-Fiber(fibercounter).PhotonPosition(i)+1) = Fiber(fibercounter).XYZ(f,1);
                YS1(f-Fiber(fibercounter).PhotonPosition(i)+1) = Fiber(fibercounter).XYZ(f,2);
                ZS1(f-Fiber(fibercounter).PhotonPosition(i)+1) = Fiber(fibercounter).XYZ(f,3);
            end
            % Obtener la máxima profundidad alcanzada por el fotón
            Fiber(fibercounter).Photondepth(i)=max(ZS1);
            % Obtener el máximo punto alcanzado sobre el eje y
            Fiber(fibercounter).PhotonYmax(i)=max(abs(YS1));
            % Obtener el máximo punto alcanzado sobre el eje x ()
            Fiber(fibercounter).PhotonXmax(i)=max(abs(XS1));
            
            % Salvar el camino del fotón en la matriz
            %[Fiber(fibercounter).Matrix] = SavePhoton(XS1,YS1,ZS1,MaxDistance,Voxel_Size(fibercounter), Number_Voxels, 1, Fiber(fibercounter).Matrix,Fiber_Separation(fibercounter),max(Fiber_Separation));            
            % Pathlength del fotón se inicia en 0
            PathLength(i)=0;
            %Se calcula el pathlength del fotón
            for n=1:size(XS1,2)-1
                %Sumando las distancias entre cada uno de los puntos de las
                %interacciones
                PathLength(i) = PathLength(i)+sqrt((XS1(n+1)-XS1(n))^2 + (YS1(n+1)-YS1(n))^2 + (ZS1(n+1)-ZS1(n))^2);
            end
            
            % Se grafica el camino del fotón
            % plot3(XS1,YS1,ZS1,'r.'), hold on
            plot3(XS1,YS1,ZS1, '.', 'color', A(k,:)), hold on
            xlabel('x [\mum]');
            ylabel('y [\mum]');
            zlabel('z [\mum]');
            title('Simulación MonteCarlo reflectancia difusa')
            if (mod(i,Num_Sources) == 0)
                k = 1;
            else
                k = k+1;
            end
            
        end

        %Guardar otros datos relevantes del camino del fotón
        %Ángulo de esparcimiento promedio
        Fiber(fibercounter).ScatAngleAverage = mean(Fiber(fibercounter).ScatAngle)*180/pi;
        %Promedio de número de eventos de esparcimiento
        Fiber(fibercounter).ScatEventsAverage = (Fiber(fibercounter).SizeScatAngle-1)/Fiber(fibercounter).PhotonCounter;
        %Intensidad (número de fotones recibidos/lanzados)
        Fiber(fibercounter).Intensity = Fiber(fibercounter).PhotonCounter/NumPhotonsLaunched;
        %Mean Pathleght
        Fiber(fibercounter).PathLength = mean(PathLength);
    end        
end

hold off

% % adds the Matrix to obtain only one plane X and Z
% boxsum = zeros(Number_Voxels,Number_Voxels,Number_Voxels);
% for i = 1:Number_Voxels
%     boxsum(:,1,:) = Fiber.Matrix(:,i,:) + boxsum(:,1,:);
% end
% 
% figure
% plot3(XS1,YS1,ZS1)
% grid

disp('Separación fibras')
%Obtener las gráficas de densidad de los fotones en el medio
for i=1:size(Fiber_Separation,2)
    i
    %Si y solo si hay al menos un fotón en el detector
    if Fiber(i).PhotonCounter>0
        %Primera figura, obtener el corte por el centro de XY y en Z=0 de
        %la densidad de fotones
        figure
        slice(Fiber(i).Matrix,Number_Voxels/2,Number_Voxels/2,0)
        shading interp
        
        %La segunda figura muestra el corte 
        %Iniciar Boxsum a ceros
        Fiber(i).Boxsum = zeros(Number_Voxels,Number_Voxels,Number_Voxels);
        %Rellenar Boxsum de 1 al número total de Voxels
        disp('# de Voxel')
        for n = 1:Number_Voxels
            n
            %Se hace el corte sobre el eje y, los valores de Matrix se
            %recorren en todo X & Z, moviendo solo Y, se hace una suma
            %acumulativa
            Fiber(i).Boxsum(:,1,:) = Fiber(i).Matrix(:,n,:)+Fiber(i).Boxsum(:,1,:);
        end
        %Se hace una copia del Boxsum
        Fiber(i).Boxsum2 = Fiber(i).Boxsum;
        %Encontrar donde Boxsum2 se hace 0
        L = find(Fiber(i).Boxsum2==0);
        %Hallar el máximo de voxel_X=50, Y & Z
        MaxVal = max(max(max(Fiber(i).Boxsum2(50,:,:))));
        %El valor del Boxsum2 en L será -MaxVal/2
        Fiber(i).Boxsum2(L)=-MaxVal/2;
        %Graficar la figura
        figure
        slice(Fiber(i).Boxsum2,1,0,0)
        shading interp
        caxis([-MaxVal/2 MaxVal]);

    end
end

%Otras gráficas relevantes
%Histograma de profunidades de penetración de los fotones
figure
subplot(2,2,1)
hist(Fiber.Photondepth,100)
xlabel('Profundidad de penetración del fotón [\mum]')
ylabel('# fotones')
title(strcat('Número de fotones detectados \newlineN=',int2str(Fiber.PhotonCounter),' fotones \newline\delta = ', num2str(ydelta*10000,'%.2f'),' \mum'),'FontWeight','bold')

% Histograma de máxima longitud alcanzada por el fotón en el eje X+
subplot(2,2,2)
hist(Fiber.PhotonXmax, 100)
xlabel('Máximo alcance del fotón eje X [\mum]')
ylabel('# fotones')
title(strcat('Número de fotones detectados \newlineN=',int2str(Fiber.PhotonCounter),' fotones'),'FontWeight','bold')

% Histograma de máxima longitud alcanzada por el fotón en el eje Y+
subplot(2,2,3)
hist(Fiber.PhotonYmax, 100)
xlabel('Máximo alcance del fotón eje Y [\mum]')
ylabel('# fotones')
title(strcat('Número de fotones detectados \newlineN=',int2str(Fiber.PhotonCounter),' fotones'),'FontWeight','bold')

% Longitud de camino medio
subplot(2,2,4)
hist(PathLength,100)
xlabel('Longitud del camino libre del fotón [\mum]')
ylabel('# fotones')
title(strcat('Número de fotones detectados \newlineN=',int2str(Fiber.PhotonCounter),' fotones \newlineCamino Libre Medio = ', num2str(Fiber.PathLength,'%.2f'),' \mum'),'FontWeight','bold')

%Mostrar las propiedades ópticas del medio y características de la fuente y
%detector
disp('<<SIMULACIÓN REFLECTANCIA DIFUSA POR MC>>')
disp(sprintf('mua       = %0.4f cm^-1', ua))
disp(sprintf('mus       = %0.4f cm^-1', us))
disp(sprintf('n_medio   = %0.4f ',n_medium))
if HG == YES %En caso de que sea x aproximaciï¿½n de HG
    disp('<-Aproximación por Henyey-Greenstein->')
    disp(sprintf('g         = %0.4f ', g))
else
    disp('<-Aproximación por Teoría de Mie->')
    disp(sprintf('n_esfera  = %0.4f ', n_sphere))
    disp(sprintf('r_esfera  = %0.4f ', sphere_radius))
    disp(sprintf('long. onda= %0.4f ', wavel))
    disp(sprintf('Ang. Mie  = %0.4f ', MieAng))
    disp(sprintf('Fase Mie  = %0.4f ', MiePhase))
    disp(sprintf('g         = %0.4f ', g))
end
disp(sprintf('Prof. Pen. Optica   = %0.4f um\n',ydelta*10000))

disp('<<FUENTE-DETECTOR>>')
disp(sprintf('n_fibra               = %0.4f ',n_fiber))
disp(sprintf('Radio_fuente          = %0.4f um',Source_Diameter))
disp(sprintf('Inclinación_fuente    = %0.4f° ',Source_Tilt_Angle))
disp(sprintf('Radio_detector        = %0.4f um',Detector_Diameter))
disp(sprintf('Inclinación_detector  = %0.4f° ',Detector_Tilt_Angle))
disp(sprintf('Separación            = %0.4f um\n',Fiber_Separation))

disp('<<DATOS SIMULACIÓN>>')
disp(sprintf('Fotones Detectados            = %0.0f ',Fiber.PhotonCounter))
disp(sprintf('Fotones Lanzados              = %0.0f ',NumPhotonsLaunched))
disp(sprintf('Intensidad                    = %0.4e ',Fiber.Intensity))
disp(sprintf('Angulo esparcimiento promedio = %0.4f° ',Fiber.ScatAngleAverage))
disp(sprintf('# eventos de esparcimiento    = %0.4f  ',Fiber.ScatEventsAverage))
disp(sprintf('Camino libre medio            = %0.4f um',Fiber.PathLength))
disp(sprintf('Energía fotón                 = %0.4e ', (hplanck*cmedium)/(lambda*10^-9) ))
disp(sprintf('Energía total                 = %0.4e ', NumPhotonsLaunched*(hplanck*cmedium)/(lambda*10^-9) ))

% figure
% slice(boxsum,1,0,1)
% shading interp    
% 
% %Mostrar la gráfica de XZ & YZ de densidad de fotones
% clear L
% Matrix2 = Fiber.Matrix;
% L = find(Matrix2==0);
% MaxVal = max(max(max(Matrix2(50,:,:))));
% Matrix2(L)=-MaxVal/2;
% figure
% slice(Matrix2,Number_Voxels/2+0.5,Number_Voxels/2+0.5,1)
% shading interp
% caxis([-MaxVal/2 MaxVal])
% 
% clear L
% boxsum2 = boxsum;
% L = find(boxsum2==0);
% MaxVal = max(max(max(boxsum2(50,:,:))));
% boxsum2(L)=-MaxVal/2;
% figure
% slice(boxsum2,1,0,1)
% shading interp
% caxis([-MaxVal/2 MaxVal])
% 
% for i = 1:Number_Voxels
%     Depth(i) = (boxsum(50,1,i));
% end
% figure
% plot(2*(MaxDistance / Number_Voxels)*(0:Number_Voxels-1),Depth)