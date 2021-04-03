function fit_scircle(Lon,Lat,HA,RadiusGuess,Units)
% 
% Example: celestial.coo.fit_scircle


RAD = 180./pi;


if nargin<5
    Units = 'deg';
    if nargin<4
        RadiusGuess = 1;
    end
end

Factor = convert.angular(Units,'rad');
Lon    = Factor.*Lon;
Lat    = Factor.*Lat;
HA     = Factor.*HA;
RadiusGuess = Factor.*RadiusGuess;


[X,Fval] = Util.fit.fminsearch_my({@radius_rms,Lon,Lat,HA},[240, 89.5, 1 0]./RAD); 

%[X,Fval] = Util.fit.fminsearch_my({@radius_rms,Lon,Lat,HA},RadiusGuess); 

Fval.*RAD
X.*RAD

-0.19

end


function [MeanDist,MeanLon,MeanLat]=radius_rms(Par,Lon,Lat,HA)
%

RAD = 180./pi;

LonC   = Par(1);
LatC   = Par(2);
Radius = Par(3);
Offset = Par(4);

[LatCalc,LonCalc]     = reckon(LatC,LonC,Radius,HA+Offset,'radians');
D = celestial.coo.sphere_dist_fast(LonCalc,LatCalc,Lon,Lat);
MeanDist = mean(D);


end



% function [MeanDist,MeanLon,MeanLat]=radius_rms(Radius,Lon,Lat,HA)
% %
% 
% [LatC,LonC]     = reckon(Lat,Lon,Radius,HA,'radians');
% [CD1,CD2,CD3] = celestial.coo.coo2cosined(LonC,LatC);
% MeanCD1 = mean(CD1);
% MeanCD2 = mean(CD2);
% MeanCD3 = mean(CD3);
% [MeanLon,MeanLat] = celestial.coo.cosined2coo(MeanCD1,MeanCD2,MeanCD3);
% 
% D = celestial.coo.sphere_dist_fast(MeanLon,MeanLat,LonC,LatC);
% MeanDist = mean(D);
% 
% end

