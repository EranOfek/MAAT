function HA=alt2ha(Alt,Dec,Phi,Units)
% Convert altitude and declnation to hour angle
% Package: celestial.coo
% Description: Given an object altitude and declination and the observer
%              latitude, return the corresponding Hour Angle.
% Input  : - Altitude [radians].
%          - Declination [radians].
%          - Geocentric latitude [radians].
%          - Input/output units 'rad' or 'deg'. Default is 'rad'.
%            Outpuoutputt is always in radians.
% Output : - Hour Angle [radians].
%            Note that there are two solutions at +/-HA, but only
%            the positive HA is returned.
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: HA=celestial.coo.alt2ha(1,1,1)
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==3)
    Units = 'rad';
end

Alt = convert.angular(Units,'rad',Alt);
Dec = convert.angular(Units,'rad',Dec);
Phi = convert.angular(Units,'rad',Phi);

% solve:
% sin(Alt)= sin(dec)*sin(phi) + cos(dec)*cos(phi)*cos(HA)
HA = abs(acos((sin(Alt) - sin(Dec).*sin(Phi))./(cos(Dec).*cos(Phi))));

HA = convert.angular('rad',Units,HA);


