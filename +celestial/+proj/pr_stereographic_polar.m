function [X,Y]=pr_stereographic_polar(Az,ZenithDist)
% Project coordinates using the Stereographic polar projection
% Package: celestial.proj
% Description: Project coordinates (longitude and latitude) using the
%              Stereographic polar projection.
%              This projection preservs angles.
% Input  : - Vector of Azimuth, in radians.
%          - Vector of Zenith-distance, in radians.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                     Nov 2004  
%    URL : http://weizmnann.ac.il/home/eofek/matlab/
% Example: [X,Y]=celestial.proj.pr_stereographic_polar(1,1)
% Reliable: 2
%--------------------------------------------------------------------------
if (nargin==2)
   % no default
elseif (nargin==3)
   error('Illigal number of argument');
end

X     = cos(Az).*tan(0.5.*ZenithDist);
Y     = sin(Az).*tan(0.5.*ZenithDist);


