function Rand=randinpolygon(Polygon,Nr,Geom)
% Generate random positions inside a polygon
% Package: Util.stat
% Description: Generate random positions inside a polygon, defined on a
%              plane or a sphere.
% Input  : - Polygon boundries [X, Y] (in case of spherical geometry
%            [long, lat] in radians.
%            Note you must connect also the first and last point!
%          - Number of random points to generate.
%          - Geometry:
%            'plane'    - polygon on a plane (default).
%            'sphere'   - polygin on a sphere.
%                         In that case the [X,Y] coordinates
%                         must be in radians.
% Output : - Coordinates [X,Y] of random points inside polygon.
% Tested : Matlab 6.5
%     By : Eran O. Ofek                    Aug 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

ColX  = 1;
ColY  = 2;  

if (nargin==2),
   Geom   = 'plane';
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

N    = 0;
Rand = zeros(0,2);
switch Geom
 case 'plane'
    %--- plane geometry ---
    Iter = 0;
    while (N<Nr),
       Iter    = Iter + 1;
       RandXY  = rand(Nr,2);
       MinX    = min(Polygon(:,ColX));
       MaxX    = max(Polygon(:,ColX));
       MinY    = min(Polygon(:,ColY));
       MaxY    = max(Polygon(:,ColY));
       RandXY(:,ColX) = RandXY(:,ColX).*(MaxX-MinX)+MinX;  
       RandXY(:,ColY) = RandXY(:,ColY).*(MaxY-MinY)+MinY;  
       In      = inpolygon(RandXY(:,ColX),RandXY(:,ColY),Polygon(:,ColX),Polygon(:,ColY));
       InInd   = find(In==1); 
       Rand    = [Rand; RandXY(InInd,:)];       
       N       = size(Rand,1);
    end
    Rand = Rand(1:Nr,:);   
 case 'sphere'
    %--- spherical geometry ---
    error('sphere geometry not available yet')
 otherwise
    error('Unknown Geom option');
end
