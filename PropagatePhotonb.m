%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funci�n para propagaci�n del foton
% Recibe:
% coeficientes �pticos (mus, mua, g)
% El paso
% Valor HG de Henyey Greenstein
% Los coeficientes de la ecuaci�n de la fuente (a, b) y su centro
% [multifuente]
% Los coeficientes de la ecuaci�n del detector (a, b)
% Posici�n de las fibras en el eje Z
% Separaci�n de las fibras
% NA
% �ngulo de inclinacion de la fuente
% �ngulo de inclinacion del detector
% indice de refracci�n de la fibra
% indice de refracci�n del medio
% M�xima Distancia entre fuente-detector
% Angulo y Fase para el c�lculo Mie
% 
% Devuelve:
% Posicion XYZ
% Contador de fotones lanzados
% �ngulo de esparcimiento
% Posicion del fot�n en la fibra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XPosition,YPosition,ZPosition,PhotonsLaunchedCounter, ScatAngle, fiber] = PropagatePhotonb(us,ua,Fix_Step,HG, g,aSource,bSource,cSource,aDetector,bDetector,Zfibers, Fiber_Separation, NA, Source_Tilt_Angle, Detector_Tilt_Angle, n_fiber, n_medium, MaxDistance, MieAng, MiePhase);

%Definici�n de variables Si/No
YES = 1;
NO = 0;

First_Photon = YES;
Photon_Detected = NO;
Photon_Alive = YES;

PhotonsLaunchedCounter = 0; % contador para el n�mero de fotones que han diso lanzados
PositionCounter = 1; % contador para la localizaci�n del fot�n
ScatAngleCounter = 1;


while Photon_Detected == NO
   
    % si es el primer fot�n que es lanzado
    if First_Photon == YES
        
        % Obtiene la coordenada XYZ de salida y la salva
        [X,Y] = StartingCoordinatesb(aSource,bSource,cSource,max(Fiber_Separation));
        Z = Zfibers;
        XPosition(PositionCounter) = X; YPosition(PositionCounter) = Y; ZPosition(PositionCounter) = Z;
        PositionCounter = PositionCounter + 1;
        
        % Obtiene la direccion
        [ux,uy,uz] = StartingDirection (NA, Source_Tilt_Angle, n_fiber, n_medium);

        % Cambia los valores iniciales
        First_Photon = NO;
        PhotonsLaunchedCounter = PhotonsLaunchedCounter + 1;
    
    % Si el fot�n no est� siendo lanzado
    else
        
        % Si se usa la aproximaci�n HG
        if HG == YES
            [ux,uy,uz,tita] = HGaprox (ux,uy,uz,g);
            
        % Si se usa la aproximacion de Mie
        else
            [ux,uy,uz,tita] = MIEaprox (ux,uy,uz,MieAng,MiePhase);
        end
        
        % Graba el �ngulo de esparcimiento
        ScatAngle(ScatAngleCounter) = tita;
        ScatAngleCounter = ScatAngleCounter + 1;
    end
    
    % Da un paso y mueve el fot�n al pr�ximo XYZ
    [X,Y,Z] = MovePhoton(XPosition(PositionCounter-1), YPosition(PositionCounter-1), ZPosition(PositionCounter-1), ux, uy, uz, Fix_Step, us, ua, Zfibers);
    %   Salva XYZ
    XPosition(PositionCounter) = X; YPosition(PositionCounter) = Y; ZPosition(PositionCounter) = Z;
    PositionCounter = PositionCounter + 1;
    
    % Si el fot�n ha ido muy lejos indica que el fot�n ya no esta "vivo" -> Photon_Alive = NO 
    if abs(XPosition(PositionCounter-1))>=MaxDistance || abs(YPosition(PositionCounter-1))>=MaxDistance || abs(ZPosition(PositionCounter-1))>=2*MaxDistance
        Photon_Alive = NO;
        
    elseif Z<=Zfibers
       
        % Ver si el foton est� bajo una de las fibras
        for fiber = 1:size(Fiber_Separation,2)
           
            if ((XPosition(PositionCounter-1) + (max(Fiber_Separation)/2 - Fiber_Separation(fiber)))^2/aDetector^2 + (YPosition(PositionCounter-1))^2/bDetector^2) <= 1      
                % Ver si el fot�n esta dentro de NA o la fibra
                
                % En caso de que sea la misma fibra fuente
                if Fiber_Separation(fiber) == 0
                    AngleOK = DetectDirection (NA, Source_Tilt_Angle, n_fiber, n_medium,ux,uy,uz);
                else
                    AngleOK = DetectDirection (NA, Detector_Tilt_Angle, n_fiber, n_medium,ux,uy,uz);
                end
                
                if AngleOK == YES
                
                    Photon_Detected=YES;
                    return    
                end
            end
        end
        
        % Si el foton no est� bajo la misma fibra o tiene un NA erroneo, 
        % entonces probar para Reflexi�n
        Reflectance = Reflected (uz, n_medium);
        if Reflectance == YES
            
            % Mover el fot�n al pr�ximo XYZ y cambiar la direcci�n Z
            uz = -uz;
            [X,Y,Z] = MovePhoton(XPosition(PositionCounter-1), YPosition(PositionCounter-1), ZPosition(PositionCounter-1), ux, uy, uz, Fix_Step, us, ua, Zfibers);
            
            % Salvar XYZ
            XPosition(PositionCounter) = X; YPosition(PositionCounter) = Y; ZPosition(PositionCounter) = Z;
            PositionCounter = PositionCounter + 1;
        else
            Photon_Alive = NO;
        end 
        
    end   
    
    
    % Photon_Detected = YES;

    % Si el fot�n est� "muerto", resetear valores
    if Photon_Alive == NO
        clear XPosition; clear YPosition; clear ZPosition;
        clear ScatAngle;
        ScatAngleCounter = 1;
        First_Photon = YES;
        PositionCounter = 1;
        Photon_Alive = YES;
    end
    
end