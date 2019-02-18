function K=gauss_grav_const(M2,M1,R,T)
% Gaussian gravitational constant for a system
% Package: celestial.Kepler
% Description: Get the analog of the Gaussian gravitational constant for
%              a system with a given primary mass, secondary mass and
%              unit distance. This program is useful in order to apply
%              kepler.m for non-solar system cases.
% Input  : - Secondary mass [solar mass units], default is 0.
%          - Primary mass [solar mass units], default is 1.
%          - Unit distance in [au], default is 1.
%          - Unit time in [day=86400 SI sec], default is 1.
% Output : - The equivalent of the Gaussian gravitational constant for
%            the chosen system. By default (K=gauss_grav_const), the
%            program return the Gaussian gravitational constant
%            (i.e., 0.017202098950000).
% See also : get_constant.m, kepler.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: K=celestial.Kepler.gauss_grav_const(0,1,1,1); % get the Gaussian gravitational constant
% Reliable: 2
%--------------------------------------------------------------------------
DefM2 = 0;
DefM1 = 1;
DefR  = 1;
DefT  = 1;
if (nargin==0)
   M2  = DefM2;
   M1  = DefM1;
   R   = DefR;
   T   = DefT;
elseif (nargin==1)
   M1  = DefM1;
   R   = DefR;
   T   = DefT;
elseif (nargin==2)
   R   = DefR;
   T   = DefT;
elseif (nargin==3)
   T   = DefT;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

K = 0.017202098950000.*R.^(3./2)./(sqrt(M1+M2).*T);
