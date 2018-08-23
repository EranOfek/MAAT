function Z=vel2shift(Vel,Units)
% Calculate the red/blue shift (z) from velocity.
% Package: AstroUtil.spec
% Description: Calculate the red/blue shift (z) from velocity.
% Input  : - Matrix of Velocities (default units are km/s).
%          - String of units of velocity. Default is 'km/s'.
% Output : - Matrix of red/blue shift (Z).
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 2000   
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: shift2vel.m
% Example: vel2shift(179875.4748)  % return redshift = 1
% Reliable: 1
%--------------------------------------------------------------------------
C = constant.c('SI')./1000;    % speed of light [km/s]

if (nargin==1)
   % do donthing - already in km/s
else
  [~,Vel] = convert.units(Units,'km/s',Vel);
end

Beta  = Vel./C;

Z = sqrt((1+Beta)./(1-Beta)) - 1;

