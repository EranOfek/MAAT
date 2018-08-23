function Vel=shift2vel(Z)
%--------------------------------------------------------------------------
% shift2vel function                                             AstroSpec
% Description: Calculate the velocity from the red/blue shift (z).
% Input  : - Vector of red/blue shift (z).
% Output : - Vector of Velocities in km/sec.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2000   
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: vel2shift.m
% Example: shift2vel(1.0)
% Reliable: 1
%--------------------------------------------------------------------------
C = get_constant('c','SI')./1000;    % speed of light [km/s]

Beta = ((Z+1).^2 - 1)./(1 + (Z+1).^2);

Vel = Beta.*C;


