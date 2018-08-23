function [Dist,PA]=sphere_dist(RA_1,Dec_1,RA_2,Dec_2,Units)
% angular distance and position angle between two points on the sphere
% Package: celestial.coo
% Description: Calculate the angular distance and position angle between
%              two points on the celestial sphere.
% Input  : - Column vector of "long." for the first point
%            [rad] or [H M S] or sexagesimal string.
%            see convertdms.m for details.
%          - Column vector of "lat." for the first point
%            [rad] or [Sign D M S] or sexagesimal string.
%          - Column vector of "long." for the second point
%            [rad] or [H M S] or sexagesimal string.
%          - Column vector of "lat." for the second point
%            [rad] or [Sign D M S] or sexagesimal string.
%          - Units of input parameters {'g'|'rad'|'deg'}, default is 'g'.
%            'g' - radians or sexagesimal.
%            'rad' - radians only.
%            'deg' - degrees only.
%            Output is always in radians.
% Output : - Vector of distances between points [radian].
%          - Vector of position angle between points [radian].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Feb 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: sphere_dist_fast.m, sphere_dist_cosd.m (and built in distance.m).
% Example: [D,P]=celestial.coo.sphere_dist([1.1;1],[0.2;0.1],[10 0 0],[1 40 0 0]);
% Reliable: 1
%------------------------------------------------------------------------------
RADIAN = 180./pi;


if (nargin==4),
   Units   = 'g';
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

switch lower(Units)
 case 'rad'
    % do nothing
 case 'g'
    RA_1  = celestial.coo.convertdms(RA_1,'gH','r');
    Dec_1 = celestial.coo.convertdms(Dec_1,'gD','r');
    RA_2  = celestial.coo.convertdms(RA_2,'gH','r');
    Dec_2 = celestial.coo.convertdms(Dec_2,'gD','r');
 case 'deg'
    RA_1   = RA_1./RADIAN;
    Dec_1  = Dec_1./RADIAN;
    RA_2   = RA_2./RADIAN;
    Dec_2  = Dec_2./RADIAN;
 otherwise
    error('Unknown Units option');
end

if (length(Dec_1)==length(Dec_2))
   % do nothing
else
   if (length(Dec_1)==1)
      Dec_1 = Dec_1.*ones(length(Dec_2),1);
   elseif (length(Dec_2)==1)
      Dec_2 = Dec_2.*ones(length(Dec_1),1);
   else
      error('Illegal vectors length');
   end
end


dRA = RA_2 - RA_1;

Dist = acos(sin(Dec_1).*sin(Dec_2) + cos(Dec_1).*cos(Dec_2).*cos(dRA));

Dist = abs(Dist);

% cosine direction method is slower
%CosD1 = cosined([RA_1, Dec_1]);
%CosD2 = cosined([RA_2, Dec_2]);
%Dist  = acos(dot(CosD1,CosD2,2));


% Checking the half sine formula: - I get the same results...
%CosDist = sin(Dec_1).*sin(Dec_2) + cos(Dec_1).*cos(Dec_2).*cos(dRA);
%%SinHalfDist = sqrt(0.5.*(1 - CosDist));
%Dist = asin(sqrt(0.5.*(1 - CosDist))).*2;



if (nargout>1),
   SinPA = sin(dRA).*cos(Dec_2)./sin(Dist);
   CosPA = (sin(Dec_2).*cos(Dec_1) - cos(Dec_2).*sin(Dec_1).*cos(dRA))./sin(Dist);
   %if (isreal(SinPA)==0 | isreal(CosPA)==0),
   %    save sd.mat RA_1, RA_2, Dec_1, Dec_2
   %end

   PA    = atan2(SinPA,CosPA);
   PA(PA<0) = 2.*pi + PA(PA<0);

   %I     = find(PA<0);
   %PA(I) = 2.*pi + PA(I);

end



