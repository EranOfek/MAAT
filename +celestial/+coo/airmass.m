function [AirMass,AzAlt,HA]=airmass(JD,RA,Dec,ObsCoo)
% Airmass from time and object and observer position
% Package: celestial.coo
% Description: Given the JD, object celestial coordinates, and observer
%              Geodetic coordinates, calculating the airmass of the
%              object.
% Input  : - Matrix in which the first column is JD (the other columns
%            are ignored).
%          - Object RA in radians or sexagesimal coordinates
%            (see convertdms.m for possible options)
%          - Object Dec in radians or sexagesimal coordinates
%            (see convertdms.m for possible options)
%          - Observer coordinates: [East Long, North Lat] in radians.
% Output : - Vector of airmass.
%          - [AZ, Alt] of object in radians.
% See Also : hardie.m; hardie_inv.m; horiz_coo.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Apr 2011
%    URL : http://wiezmann.ac.il/home/eofek/matlab/
% Example: [AirMass,AzAlt]=celestial.coo.airmass(2451545+[0;1],1,1,[1 1]);
% Reliable: 2
%--------------------------------------------------------------------------

RA  = celestial.coo.convertdms(RA,'gH','r');
Dec = celestial.coo.convertdms(Dec,'gD','R');

AzAlt = celestial.coo.horiz_coo([RA, Dec],JD, ObsCoo,'h');   
AirMass = celestial.coo.hardie(pi./2 - AzAlt(:,2));

if (nargout>2)
    LST = celestial.time.lst(JD,ObsCoo(:,1));
    HA  = LST.*2.*pi - RA;
end

