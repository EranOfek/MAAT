function [DomeAz,Az,Alt]=dome_az(HA,Dec,ObsLat,DomeRad,MountRad)
%--------------------------------------------------------------------------
% dome_az function                                                   ephem
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    May 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DomeAz,Az,Alt]=dome_az(0,0,32./RAD,5,1)
% Reliable: 
%--------------------------------------------------------------------------

need to take z into account ...

SinAlt = sin(Dec).*sin(ObsLat) + cos(Dec).*cos(HA).*cos(ObsLat);
CosAlt = sqrt(1-SinAlt.*SinAlt);

SinAz  = (-cos(Dec).*sin(HA))./CosAlt;
CosAz  = (sin(Dec).*cos(ObsLat) - cos(Dec).*cos(HA).*sin(ObsLat))./CosAlt;
Az     = atan2(SinAz,CosAz);
Alt    = asin(SinAlt);

% X/Y position of telescope-delination axis
Xt = MountRad.*cos(ObsLat).*sin(HA);
Yt = MountRad.*sin(ObsLat).*cos(HA);

Xaz = cos(Az);
Yaz = sin(Az);

Xdome = Xaz-Xt;
Ydome = Yaz-Yt;

DomeAz = atan2(Ydome,Xdome);

