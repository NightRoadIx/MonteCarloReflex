%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Reflectance] = Reflected(uz, n);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Reflectance] = Reflected(uz, n)

YES = 1;
NO = 0;

IncidentAngle = acos(abs(uz));
TransmitAngle = asin(sin(IncidentAngle)*n/1);
Reflectance = 1/2*((sin(IncidentAngle-TransmitAngle)/sin(IncidentAngle+TransmitAngle))^2 + (tan(IncidentAngle-TransmitAngle)/tan(IncidentAngle+TransmitAngle))^2);
if rand <= Reflectance
    Reflectance = YES;
else
    Reflectance = NO;
end