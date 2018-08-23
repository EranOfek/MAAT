function E=true2eccentric_anomaly(Nu,ecc)
% True anomaly to eccentric anomaly
% Package: celestial.Kepler
% Description: Convert true anomaly to eccentric anomaly for elliptic
%              orbit.
% Input  : - True anomaly
%          - Eccentricity.
% Output : - Eccentric anomaly.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: E=celestial.Kepler.true2eccentric_anomaly(1,0.1);
% Reliable: 2
%--------------------------------------------------------------------------

E = 2.*atan(sqrt((1-ecc)./(1+ecc)).*tan(0.5.*Nu));
