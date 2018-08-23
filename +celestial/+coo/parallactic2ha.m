function [HA1,HA2]=parallactic2ha(PA,Dec,Lat)
% Convert parallactic angle and declinatio to hour angle
% Package: celestial.coo
% Description: Convert parallactic angle, declination and latitude to
%              hour angle. Note that there are two solutions, and the
%              function will return both.
% Input  : - Parallactic angle [rad].
%          - Declination [rad].
%          - Observer latitude [rad].
% Output : - Hour angle first solution [rad].
%          - Hour angle second solution [rad].
% Tested : Matlab R 2011b
%     By : Eran O. Ofek                    Nov 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [HA1,HA2]=celestial.coo.parallactic2ha(1,1,1);
% Reliable: 2
%--------------------------------------------------------------------------

HA2  = +2.*atan((tan(Lat./2).^2.*tan(PA./2).^2 - tan(Dec./2).^2.*tan(PA./2).^2 - tan(Dec./2).^2.*tan(Lat./2).^2 + (-(2.*tan(PA./2) + 2.*tan(Lat./2).^2.*tan(PA./2) - tan(Dec./2).^2.*tan(Lat./2).^2 + tan(Dec./2).^2.*tan(PA./2).^2 - tan(Lat./2).^2.*tan(PA./2).^2 + tan(Dec./2).^2 - tan(Lat./2).^2 + tan(PA./2).^2 - 2.*tan(Dec./2).^2.*tan(PA./2) - 2.*tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2) - tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2).^2 + 1).*(2.*tan(PA./2) + 2.*tan(Lat./2).^2.*tan(PA./2) + tan(Dec./2).^2.*tan(Lat./2).^2 - tan(Dec./2).^2.*tan(PA./2).^2 + tan(Lat./2).^2.*tan(PA./2).^2 - tan(Dec./2).^2 + tan(Lat./2).^2 - tan(PA./2).^2 - 2.*tan(Dec./2).^2.*tan(PA./2) - 2.*tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2) + tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2).^2 - 1)).^(1./2) + tan(Dec./2).^2 - tan(Lat./2).^2 - tan(PA./2).^2 + tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2).^2 + 1)./(4.*(tan(Dec./2).*tan(PA./2) + tan(Lat./2).*tan(PA./2) - tan(Dec./2).*tan(Lat./2).^2.*tan(PA./2) - tan(Dec./2).^2.*tan(Lat./2).*tan(PA./2))));
HA1  = -2.*atan((tan(Dec./2).^2.*tan(Lat./2).^2 + tan(Dec./2).^2.*tan(PA./2).^2 - tan(Lat./2).^2.*tan(PA./2).^2 + (-(2.*tan(PA./2) + 2.*tan(Lat./2).^2.*tan(PA./2) - tan(Dec./2).^2.*tan(Lat./2).^2 + tan(Dec./2).^2.*tan(PA./2).^2 - tan(Lat./2).^2.*tan(PA./2).^2 + tan(Dec./2).^2 - tan(Lat./2).^2 + tan(PA./2).^2 - 2.*tan(Dec./2).^2.*tan(PA./2) - 2.*tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2) - tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2).^2 + 1).*(2.*tan(PA./2) + 2.*tan(Lat./2).^2.*tan(PA./2) + tan(Dec./2).^2.*tan(Lat./2).^2 - tan(Dec./2).^2.*tan(PA./2).^2 + tan(Lat./2).^2.*tan(PA./2).^2 - tan(Dec./2).^2 + tan(Lat./2).^2 - tan(PA./2).^2 - 2.*tan(Dec./2).^2.*tan(PA./2) - 2.*tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2) + tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2).^2 - 1)).^(1./2) - tan(Dec./2).^2 + tan(Lat./2).^2 + tan(PA./2).^2 - tan(Dec./2).^2.*tan(Lat./2).^2.*tan(PA./2).^2 - 1)./(4.*(tan(Dec./2).*tan(PA./2) + tan(Lat./2).*tan(PA./2) - tan(Dec./2).*tan(Lat./2).^2.*tan(PA./2) - tan(Dec./2).^2.*tan(Lat./2).*tan(PA./2))));

