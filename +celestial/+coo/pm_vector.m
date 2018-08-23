function M=pm_vector(RA,Dec,MuRA,MuDec,V,Par)
%--------------------------------------------------------------------------
% pm_vector function                                                 ephem
% Description: Return the space motion vector given proper motion, parallax
%              and radial velocity.  
%              Obsolete: use pm2space_motion instead.
% Input  : - R.A. in radians.
%          - Dec. in radians.
%          - Proper motion in R.A. in radians/century.
%          - Proper motion in Dec. in radians/century.
%          - Radial velocity (V) in au/century (1km/sec = 21.095 au/century),
%            measured positively away frm the Earth.
%            Default is 0.
%          - Parallax in radians. Default is 0.
% Output : - Space motion vector (m) of the star expressed
%            in radians per cenury.
%            Referred to the standard equator and the coordinates equinox.
%            If the input parameters are vectors, then the result
%            is a matrix: each column for each vector element.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: M=pm_vector(1,1,1e-9,1e-9,10,1e-9);
% Reliable: 2
%--------------------------------------------------------------------------
if (nargin==4)
   V   = zeros(size(RA));
   Par = zeros(size(RA));
elseif (nargin==5)
   Par = zeros(size(RA));
elseif (nargin==6)
   % do nothing
else
   error('Illegal number of input arguments');
end

Mx = -MuRA.*cos(Dec).*sin(RA) - MuDec.*sin(Dec).*cos(RA) + V.*Par.*cos(Dec).*cos(RA);
My =  MuRA.*cos(Dec).*cos(RA) - MuDec.*sin(Dec).*sin(RA) + V.*Par.*cos(Dec).*sin(RA);
Mz =                            MuDec.*cos(Dec)          + V.*Par.*sin(Dec);

SizeMx = size(Mx);
if (SizeMx(1)>=SizeMx(2)),
   Mx = Mx.';
   My = My.';
   Mz = Mz.';
end

M = [Mx;My;Mz];

   
