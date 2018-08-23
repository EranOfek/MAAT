function planar_sundial(Lat,GnomDec,ZenithDist,Length)
% Calculate and plot a planar sundial
% Package: celestial.SolarSys
% Description: Calculate and plot a planar sundial.
% Input  : - Latitude of sundial [radians].
%          - The azimuth of the perpendicular to the sundial
%            plane, measured from the north, eastward [radians].
%          - The zenith distance of the direction defined by
%            the sundial straight stylus [radians].
%          - The length of the straight stylus [length units].
% Output : null
% Plot   : plot a sundial
%          The x-axis is the horizonatal axis;
%          the y-axis coincides with the line of greatest slope of
%          the sundial plane; and the red + mark the position of
%          the stylus footprint.
%          The line at the top mark the stylus length.
% Reference: Astronomical Algo. Meeus.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Aug 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.SolarSys.planar_sundial(1,1,1,1)
% Reliable: 2
%------------------------------------------------------------------------------

RAD = 180./pi;

MinAlt = 5./RAD;    % minimum altitude above the horizon to plot
MinHA  = -8;
MaxHA  = 8;
MinDec = -23.5;
MaxDec = 23.5;

GnomDec = GnomDec - pi;

TimeZoneOffset = +20./60;   % hours (different between Long and time zone)

HA      = (TimeZoneOffset+(MinHA:0.1:MaxHA).')./24 .*2.*pi;
HAo     = (TimeZoneOffset+(MinHA:1.0:MaxHA).')./24 .*2.*pi;
Dec     = [MinDec; 0; MaxDec]./RAD;
DecGrid = (MinDec:0.5:MaxDec).'./RAD;

figure;
hold on;
for Id=1:1:length(Dec),
   SinAlt = sin(Dec(Id)).*sin(Lat) + cos(Dec(Id)).*cos(HA).*cos(Lat);
   Alt    = asin(SinAlt);
   
   [P,Q,Nx,Ny,X,Y]=sundial_pars(Lat,GnomDec,ZenithDist,Length,HA,Dec(Id));
  
   Iab    = find(Q>MinAlt & Alt>MinAlt);

   plot(X(Iab),Y(Iab),'-','LineWidth',1,'Color',[0.5 0.5 0.5]);
end

for Ih=1:1:length(HAo),
   SinAlt = sin(DecGrid).*sin(Lat) + cos(DecGrid).*cos(HAo(Ih)).*cos(Lat);
   Alt    = asin(SinAlt);

   [P,Q,Nx,Ny,X,Y]=sundial_pars(Lat,GnomDec,ZenithDist,Length,HAo(Ih),DecGrid);

   Iab    = find(Q>MinAlt & Alt>MinAlt);

   plot(X(Iab),Y(Iab),'k-','LineWidth',2);
end

plot(0,0,'r+','MarkerSize',10);
axis([-4 4 -4 4]);
H=plotline(-0.5.*Length,3.8,Length,0);
axis square;
axis off;



%-----------------------------------------------------------
function [P,Q,Nx,Ny,X,Y]=sundial_pars(Lat,GnomDec,ZenithDist,Length,HA,Dec)
%-----------------------------------------------------------
P = sin(Lat).*cos(ZenithDist) - cos(Lat).*sin(ZenithDist).*cos(GnomDec);
Q = sin(GnomDec).*sin(ZenithDist).*sin(HA) + ...
    (cos(Lat).*cos(ZenithDist) + ...
     sin(Lat).*sin(ZenithDist).*cos(GnomDec)).*cos(HA) + ...
    P.*tan(Dec);
Nx= cos(GnomDec).*sin(HA) - ...
    sin(GnomDec).*(sin(Lat).*cos(HA) - cos(Lat).*tan(Dec));
Ny= cos(ZenithDist).*sin(GnomDec).*sin(HA) - ...
    (cos(Lat).*sin(ZenithDist) - ...
     sin(Lat).*cos(ZenithDist).*cos(GnomDec)).*cos(HA) - ...
    (sin(Lat).*sin(ZenithDist) + ...
     cos(Lat).*cos(ZenithDist).*cos(GnomDec)).*tan(Dec);

X = Length.*Nx./Q;
Y = Length.*Ny./Q;




   
