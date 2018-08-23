function Nu=eccentric2true_anomaly(E,ecc)
% Convert Eccentric anomaly to true anomaly
% Package: celestial.Kepler
% Description: Convert Eccentric anomaly to true anomaly.
% Input  : - Eccentric anomaly [radians].
%          - Eccentricity.
% Output : - True anomaly [radians].
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Nu=celestial.Kepler.eccentric2true_anomaly(1,0.1)
% Reliable: 2
%--------------------------------------------------------------------------

Nu = 2.*atan(sqrt((1+ecc)./(1-ecc)).*tan(0.5.*E));

