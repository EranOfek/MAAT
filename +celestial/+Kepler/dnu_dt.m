function [dnudt,drdt]=dnu_dt(n,e,a,E)
% Calculate dnu/dt and dr/dt for elliptical orbit
% Package: celestial.Kepler
% Description: Calculate dnu/dt and dr/dt for elliptical orbit using
%              the Kepler Equation.
% Input  : - Mean motion [rad/day].
%          - Eccentricity.
%          - Semi major axis [au].
%          - Eccentric anomaly [rad].
% Output : - dnu/dt - the time derivative of the the true anomaly
%            [rad/day].
%          - dr/dt - the time derivative of the radius vector [au/day].
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [dnudt,drdt]=celestial.Kepler.dnu_dt(0.017,0.2,1,0)
% Reliable: 2
%--------------------------------------------------------------------------

dnudt = n.*sqrt(1-e.^2)./((1-e.*cos(E)).^2);
drdt  = n.*a.*e.*sin(E)./(1-e.*cos(E));
