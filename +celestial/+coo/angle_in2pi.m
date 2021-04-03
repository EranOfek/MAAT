function Long=angle_in2pi(Long,Period)
% Convert an angle to the 0 to 2*pi range
% Package: celestial.coo
% Description: Convert an angle to the range 0 to 2.*pi.
% Input  : - Matrix of angles.
%          - Period of angles. Default is 2.*pi.
% Output : - Angle in the allowed range.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.coo.angle_in2pi([-370;-350;-10;10;350;370;730],360)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
   Period = 2.*pi;
end

Long = mod(Long,Period);

%Long = (Long./Period - floor(Long./Period)).*Period;

