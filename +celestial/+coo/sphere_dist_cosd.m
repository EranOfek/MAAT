function [Dist]=sphere_dist_cosd(CD1,CD2)
% Angular distance between a set of two cosine vector directions.
% Package: celestial.coo
% Description: Calculate the angular distance between a set of two
%              cosine vector directions.
%              This should be used instead of sphere_dist_fast.m only
%              if you have the cosine vectors.
% Input  : - A 3 column matrix of cosine directions, row per vector.
%          - A 3 column matrix of cosine directions, row per vector.
% Output : - Vector of distances between points [radian].
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: sphere_dist.m, sphere_dist_fast.m (and built in distance.m).
% Example: D=sphere_dist_cosd(CD1,CD2);
% Reliable: 1
%--------------------------------------------------------------------------

Dist = acos(sum(CD1.*CD2,2)./(sqrt(sum(CD1.^2,2).*sum(CD2.^2,2)))  );
%Dist = acos(dot(CD1,CD2,2)./(sqrt(sum(CD1.^2,2).*sum(CD2.^2,2)))  );
