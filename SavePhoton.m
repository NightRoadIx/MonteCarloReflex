%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Matrix] = SavePhoton(XSave,YSave,ZSave,MaxDistance,Voxel_Size, Number_Voxels, W, Matrix)
% Esta función regresa la Matrix con una actualización de donde el fotón ha
% viajado 
%
% Argumentos:
%    - XSave,YSave,ZSave: Las coordenadas de donde ha estado el fotón
%    - MaxDistance: La más lejana distancia en la que se sigue al fotón
%    - Voxel_Size: Tamaño del voxel en um
%    - Number_Voxels: # de voxeles dentro del medio
%    - W: Peso del fotón
%    - Matrix: Matriz que sigue la posición del fotón
%
% Regresa:
%    - Matrix: La matriz que sigue la posición del fotón
%
% Desarrollado por: Roberto Reif / Boston University 
% Última actualización: Julio 13, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Matrix] = SavePhoton(XSave,YSave,ZSave,MaxDistance,Voxel_Size, Number_Voxels, W, Matrix,FSep,MaxSep)

NO = 0;
YES = 1;

MaxPoints = size(XSave,2);

% Obtiene l voxel y lo salva para la coordenada de salida de la fibra
[xvoxbm, yvoxbm, zvoxbm] = GetVoxel (XSave(1),YSave(1),ZSave(1),MaxDistance,Voxel_Size,Number_Voxels,FSep,MaxSep);
Matrix(xvoxbm,yvoxbm,zvoxbm) = Matrix(xvoxbm,yvoxbm,zvoxbm) + W;

for i = 1:MaxPoints-1
   
    % Obtiene los puntos de inicio y fin
    XStart = XSave(i); YStart = YSave(i); ZStart = ZSave(i);
    XEnd = XSave(i+1); YEnd = YSave(i+1); ZEnd = ZSave(i+1);
    
    % La distancia entre los puntos inicial y final
    maxDist = sqrt((XStart-XEnd)^2 + (YStart-YEnd)^2 + (ZStart-ZEnd)^2);
   
    % El vector unidad a la dirección entre los puntos inicial y final
    ux = (XEnd - XStart)/(maxDist+eps);
    uy = (YEnd - YStart)/(maxDist+eps);
    uz = (ZEnd - ZStart)/(maxDist+eps);
    
    % Los puntos iniciales
    xbm = XSave(i); % xbm es 'x antes de moverse' y xam es 'x después'
    ybm = YSave(i);
    zbm = ZSave(i);
    
    % Obtiene el voxel bm y el voxel final
    [xvoxbm, yvoxbm, zvoxbm] = GetVoxel (xbm,ybm,zbm,MaxDistance,Voxel_Size,Number_Voxels,FSep,MaxSep);
    [xvoxend, yvoxend, zvoxend] = GetVoxel (XEnd,YEnd,ZEnd,MaxDistance,Voxel_Size,Number_Voxels,FSep,MaxSep);
    if zvoxend<=0
        ZSave
    end
    
    
    ReachEndPoint = NO;
    
    while ReachEndPoint==NO
        % Obtiene el punto y voxel de am
        [xam, yam, zam] = MovePhoton (xbm, ybm, zbm, ux, uy, uz, YES,1/Voxel_Size*10000,0,0);
        [xvoxam, yvoxam, zvoxam] = GetVoxel (xam,yam,zam,MaxDistance,Voxel_Size,Number_Voxels,FSep,MaxSep);
        
        % si aún no ha alcanzado la máxima distancia
        if sqrt((XStart-xam)^2 + (YStart-yam)^2 + (ZStart-zam)^2)<maxDist-(1e-12)% Esta es una pequeña corrección de MovePhoton
            % y si los voxels son diferentes los salva a la matriz
            if abs(xvoxam-xvoxbm)+abs(yvoxam-yvoxbm)+abs(zvoxam-zvoxbm)~=0
                Matrix(xvoxam,yvoxam,zvoxam) = Matrix(xvoxam,yvoxam,zvoxam) + W;    
            end
            % mover am a bm
            xvoxbm = xvoxam; yvoxbm = yvoxam; zvoxbm = zvoxam;
            xbm = xam; ybm = yam; zbm = zam;
        % Si llegaron más lejos que la máxima distancia    
        else
            %y si es un voxel diferente lo salva
            if abs(xvoxbm-xvoxend)+abs(yvoxbm-yvoxend)+abs(zvoxbm-zvoxend)~=0
                Matrix(xvoxend,yvoxend,zvoxend) = Matrix(xvoxend,yvoxend,zvoxend) + W;  
            end
            ReachEndPoint = YES;
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x1,y1,z1] = GetVoxel (X,Y,Z,TrackableRadius,Voxel_Size,Number_Voxels)
% Esta función regresa los valores x1, y1, & z1 del voxel dentro de la
% Matriz. Los valores son números redondeados mayores que 1
%
% Argumentos:
%    - X,Y,Z: localización del fotón
%    - TrackableRadius: El radio después de que el fotón se pierde
%    - Voxel_Size: tamaño del voxel en um
%    - Number_Voxels: número de voxels en la Matrix
%
% Regresa:
%    - x1,y1,z1: Las coordenadas del voxel. Los valores son números
%    redondeados
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1,y1,z1] = GetVoxel (X,Y,Z,MaxDistance,Voxel_Size,Number_Voxels,FSep,MaxSep)

% Obtiene los valores del voxel en números redondeados. X & Y deben de ser
% movidos por lo que son > 0
x1 = min(ceil(((X+MaxDistance+(MaxSep-FSep)/2)/Voxel_Size) + eps),Number_Voxels);
y1 = min(ceil(((Y+MaxDistance+(MaxSep-FSep)/2)/Voxel_Size) + eps), Number_Voxels);
z1 = min(ceil((Z)/Voxel_Size+eps), Number_Voxels);