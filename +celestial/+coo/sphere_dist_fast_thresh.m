function [Dist,PA]=sphere_dist_fast_thresh(RA_1,Dec_1,RA_2,Dec_2,Thresh)
%--------------------------------------------------------------------------
% sphere_dist_fast_thresh function                                   ephem
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
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Feb 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: sphere_dist.m (and built in distance.m).
% Example: D=sphere_dist_fast_thresh(RA1,Dec1,RA2,Dec2);
% Reliable: 
%--------------------------------------------------------------------------

Dist = acos(sin(Dec_1).*sin(Dec_2) + cos(Dec_1).*cos(Dec_2).*cos(RA_1-RA_2));

if (nargout>1),
   I = find(Dist<Thresh);
   I1 = min(I,numel(RA_1));
   I2 = min(I,numel(RA_2));
   PA = zeros(size(Dist)).*NaN;
   
   dRA = RA_1(I1)-RA_2(I2);
   SinPA = sin(dRA).*cos(Dec_2(I2))./sin(Dist(I));
   CosPA = (sin(Dec_2(I2)).*cos(Dec_1(I1)) - cos(Dec_2(I2)).*sin(Dec_1(I1)).*cos(dRA))./sin(Dist(I));
  
   PA(I)    = atan2(SinPA,CosPA);
   PA(PA<0) = 2.*pi + PA(PA<0);

end
