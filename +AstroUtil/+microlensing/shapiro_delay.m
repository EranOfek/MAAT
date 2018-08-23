function Dt=shapiro_delay(Beta,M)
% Calculate the Shapiri time delay approximation (beta>ThetaE)
% Package: AstroUtil.microlensing
% Description: Calculate the Shapiri time delay approximation assuming
%              the angle between the lens and the source is much larger
%              than the Einstein radius.
% Input  : - Angle (radians) between source and lens.
%          - Lens mass [solar mass]. Default is 1.
% Output : - Time delay [seconds]
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Dt=AstroUtil.microlensing.shapiro_delay(0.001);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<2)
    M = 1;
end

Dt = -2.*constant.G.*constant.SunM.*M./(constant.c.^3) .* log(1 - cos(Beta));
