function [Dist,Ang,PA]=sphere_dist_fast(RA_1,Dec_1,RA_2,Dec_2)
%--------------------------------------------------------------------------
% sphere_dist_fast function                                          ephem
% Description: Calculate the angular distance between two points on the
%              celestial sphere. See sphere_dist.m (and built in distance.m)
%                for a more general function. This function is ~10 time
%              faster than sphere_dist.m, but it works only with radians
%                and calculate only the distance.
% Input  : - Matrix of logitudes for the first point [radian].
%          - Matrix of latitudes for the first point [radian].
%          - Matrix of logitudes for the second point [radian].
%          - Matrix of latitudes for the second point [radian].
% Output : - Matrix of distances between points [radian].
%          - Matrix of position angles between points [radian].
%            Measured westward. Take 2*pi-PA to get the PA measured
%            Eastward.
%          - Matrix of P.A. of the first point relative to the second point
%            (eastwrad from the north).
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Feb 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: sphere_dist.m, sphere_dist_cosd.m (and built in distance.m).
% Example: D=celestial.coo.sphere_dist_fast(RA1,Dec1,RA2,Dec2);
% Reliable: 1
%--------------------------------------------------------------------------

Dist = acos(sin(Dec_1).*sin(Dec_2) + cos(Dec_1).*cos(Dec_2).*cos(RA_1-RA_2));
%Dist = acos(sin(Dec_1).*sin(Dec_2) + sqrt(1-sin(Dec_1).^2).*sqrt(1-sin(Dec_2).^2).*cos(RA_1-RA_2));  % this is more accurate


if (nargout>1)
   dRA = RA_1-RA_2;
   SinPA = sin(dRA).*cos(Dec_2)./sin(Dist);
   CosPA = (sin(Dec_2).*cos(Dec_1) - cos(Dec_2).*sin(Dec_1).*cos(dRA))./sin(Dist);
  
   Ang    = atan2(real(SinPA),real(CosPA));
   %PA(PA<0) = 2.*pi + PA(PA<0);
   
   I     = find(Ang<0);
   Ang(I) = 2.*pi + Ang(I);
   
   if nargout>2
        PA    = atan2(real(SinPA),-real(CosPA));
        PA(PA<0) = PA(PA<0) + 2.*pi;
   end

end
