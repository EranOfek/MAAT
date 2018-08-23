function [X,Y]=projectcoo(Long,Lat,ProjectionType,R,Par1,Par2)
%--------------------------------------------------------------------------
% projectcoo function                                             AstroMap
% Description: Project coordinates from longitude and latitude
%              to X/Y using a specified projection.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Projection Type
%            'a' - Aitoff (default)
%            'm' - Mollweide equal-area
%            'h' - Hammer
%            'p' - Parabolic
%            's' - Sinusoidal
%            'l' - Lambert equal-area cylindrical
%            'b' - Behrmann equal-area cylindrical
%            't' - Tristan Edwards cylindrical
%            'P' - Peters cylindrical
%            'G' - Gall Orthographic cylindrical
%            'B' - Balthasart cylindrical
%            'c' - General cylindrical, opt par 1 is [Stand_Long, Stand_Lat] (radians)
%            'C' - Cassini
%            'x' - XY projection, no transformation.
%            'r' - polar projection (from north pole).
%            'A' - Albers equal-area, Par 1 is CenCoo and 2 is for ParLat (radians)
%            'g' - Gnomonic nonconformal projection. Par1 is [Long_cen, Lat_cen]
%            'M' - Mercator projection. Par1 is long_cen.
%            'o' - Bonne projection. Par1 is [central_long, standard_parallel]
%            'S' - Stereographic projection. Par1 is [central_long, standard_parallel]
%            #.# - Conic projection. where #.# is height of apex.
%          - Radius scale parameters, default is 1.
%          - Optional parameters 1.
%          - Optional parameters 2.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.2
%     By : Eran O. Ofek                    Jul 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
import celestial.proj.*

RADIAN = 180./pi;

if (nargin~=4)
   R = 1;
end

if (ProjectionType=='h')
   [X,Y]=pr_hammer(Long,Lat);
elseif (ProjectionType=='a')
   [X,Y]=pr_aitoff(Long,Lat,R);
elseif (ProjectionType=='m')
   [X,Y]=pr_mollweide(Long,Lat,R);
elseif (ProjectionType=='p')
   [X,Y]=pr_parabolic(Long,Lat,R);
elseif (ProjectionType=='s')
   [X,Y]=pr_sinusoidal(Long,Lat,R);
elseif (ProjectionType=='l')
   [X,Y]=pr_cylindrical(Long,Lat,R,[0 0]);
elseif (ProjectionType=='b')
   [X,Y]=pr_cylindrical(Long,Lat,R,[0 30/RADIAN]);
elseif (ProjectionType=='t')
   [X,Y]=pr_cylindrical(Long,Lat,R,[0 37.383/RADIAN]);
elseif (ProjectionType=='P')
   [X,Y]=pr_cylindrical(Long,Lat,R,[0 44.138/RADIAN]);
elseif (ProjectionType=='G')
   [X,Y]=pr_cylindrical(Long,Lat,R,[0 45/RADIAN]);
elseif (ProjectionType=='B')
   [X,Y]=pr_cylindrical(Long,Lat,R,[0 50/RADIAN]);
elseif (ProjectionType=='c')
   if (nargin==5)
      [X,Y]=pr_cylindrical(Long,Lat,R,Par1);
   else   
      [X,Y]=pr_cylindrical(Long,Lat,R);
   end
elseif (ProjectionType=='C')
   [X,Y]=pr_cassini(Long,Lat,R);
elseif (ProjectionType=='x')
   [X,Y]=pr_xy(Long,Lat,R);
elseif (ProjectionType=='r')
   [X,Y]=pr_polar(Long,Lat,R);
elseif (ProjectionType=='A')
   if (nargin==6)
      [X,Y]=pr_albers(Long,Lat,R,Par1,Par2);
   elseif (nargin==5)
      [X,Y]=pr_albers(Long,Lat,R,Par1);
   else
      [X,Y]=pr_albers(Long,Lat,R);
   end
elseif (ProjectionType=='g')
   if (nargin==5)
      [X,Y]=pr_gnomonic(Long,Lat,R,Par1);
   elseif (nargin==4)
      [X,Y]=pr_gnomonic(Long,Lat,R,[0 0]);
   elseif (nargin==3)
      [X,Y]=pr_gnomonic(Long,Lat,1,[0 0]);
   else
      error('Illigal number of parameters')
   end
elseif (ProjectionType=='M')
   if (nargin==5)
      [X,Y]=pr_mercator(Long,Lat,R,Par1);
   elseif (nargin==4)
      [X,Y]=pr_mercator(Long,Lat,R,0);
   elseif (nargin==3)
      [X,Y]=pr_mercator(Long,Lat,1,0);
   else
      error('Illigal number of parameters')
   end
elseif (ProjectionType=='o')
   if (nargin==5)
      [X,Y]=pr_bonne(Long,Lat,R,Par1);
   elseif (nargin==4)
      [X,Y]=pr_bonne(Long,Lat,R,[0 pi./4]);
   elseif (nargin==3)
      [X,Y]=pr_bonne(Long,Lat,1,[0 pi./4]);
   else
      error('Illigal number of parameters')
   end
elseif (ProjectionType=='S')
   if (nargin==5)
      [X,Y]=pr_stereographic(Long,Lat,R,Par1);
   elseif (nargin==4)
      [X,Y]=pr_stereographic(Long,Lat,R,[0 0]);
   elseif (nargin==3)
      [X,Y]=pr_stereographic(Long,Lat,1,[0 0]);
   else
      error('Illigal number of parameters')
   end
elseif ~ischar(ProjectionType)
   H = ProjectionType;
   if (nargin==4)
      [X,Y]=pr_conic(Long,Lat,R,H);
   elseif (nargin==3)
      [X,Y]=pr_conic(Long,Lat,1,H);
   else
      error('Illigal number of parameters')
   end
else
   error('Unknown projection type');
end


