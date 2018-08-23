function [Dec1,Dec2]=altha2dec(Alt,HA,Phi,Units)
% Convert altitude and hour angle to declination
% Package: celestial.coo
% Description: Given Altitude and Hour Angle of an object and the observer
%              latitude, calculate the object Declination.
%              There may be up to two solutions for the Declination.
% Input  : - Altitude [radians].
%          - Hour Angle [radians].
%          - Observer geocentric latitude [radians].
%          - Input units 'rad' or 'deg'. Default is 'rad'.
%            Output is always in radians.
% Output : - First Declination solution [radians].
%            NaN if no solution.
%          - Second Declination solution [radians].
%            NaN if no solution.
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Alt=40;  Phi=32; HA= (-180:1:180)';
%          [Dec1,Dec2]=celestial.coo.altha2dec(Alt,HA,Phi,'deg')
%          plot(HA,Dec1.*RAD);hold on; plot(HA,Dec2.*RAD,'r-') 
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==3),
    Units = 'rad';
end
switch lower(Units)
    case 'rad'
        % do nothing
    case 'deg'
        % convert deg to rad
        InvRAD = pi./180;
        
        Alt = Alt.*InvRAD;
        HA  = HA.*InvRAD;
        Phi = Phi.*InvRAD;
    otherwise
        error('Unknown Units option');
end

% solve:
% sin(Alt)= sin(dec)*sin(phi) + cos(dec)*cos(phi)*cos(HA)
% -> sin(Alt) - sin(dec)*sin(phi) = cos(dec)*cos(phi)*cos(HA)
% take the power of two
% and solve cubic equation in sin(dec).
CosPhiHA2 = (cos(Phi)*cos(HA)).^2;
A  = sin(Phi).^2 + CosPhiHA2;
B  = -2.*sin(Alt).*sin(Phi);
C = sin(Alt).^2 - CosPhiHA2;

Dec1 = asin((-B + sqrt(B.^2 - 4.*A.*C))./(2.*A));
Dec2 = asin((-B - sqrt(B.^2 - 4.*A.*C))./(2.*A));
Dec1(imag(Dec1)~=0)=NaN;
Dec2(imag(Dec2)~=0)=NaN;


