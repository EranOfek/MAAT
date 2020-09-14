function [OffsetLong,OffsetLat,Dist,PA]=sphere_offset(RA1,Dec1,RA2,Dec2,Units,OffsetType)
%--------------------------------------------------------------------------
% sphere_offset function                                             ephem
% Description: Calculate the offset needed to move from a point on the
%              celesial sphere to a second point on the celestial sphere,
%              along longitide (small circle) and latitude (great circle).
%              The needed offsets depends in which axis the offset is done
%              first (longitude or latitude - 'rd' and 'dr' options, 
%              respectively).
% Input  : - Column vector of "long." for the first point
%            [rad] or [H M S] or sexagesimal string.
%          - Column vector of "lat." for the first point
%            [rad] or [Sign D M S] or sexagesimal string.
%          - Column vector of "long." for the second point
%            [rad] or [H M S] or sexagesimal string.
%          - Column vector of "lat." for the second point
%            [rad] or [Sign D M S] or sexagesimal string.
%          - Units of input parameters {'rad'|'deg'}, default is 'rad'.
%            Output is always in radians.
%          - Type of offset:
%            'rd'   : Move in longitude (e.g., RA) followed by move
%                     in latitude (e.g., Dec.). Default.
%            'dr'   : First move in latitude, then in longitude.
% Output : - Vector of offsets in longitude [rad].
%          - Vector of offsets in latitude [rad].
%          - Vector of distances between the two points (along great
%            circle; see sphere_dist.m).
%          - Vector of position angles between the two points.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [OffsetLong,OffsetLat,Dist,PA]=sphere_offset([10 0 0],[1 60 0 0],[10 0 15],[1 60 2 56],[],'rd');
% Reliable: 1
%--------------------------------------------------------------------------

RAD   = 180./pi;
DefUnits      = 'rad';
DefOffsetType = 'rd';
if (nargin==4)
   Units       = DefUnits;
   OffsetType  = DefOffsetType;
elseif (nargin==5)
   OffsetType  = DefOffsetType;
elseif (nargin==6)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Units)==1)
   Units = 'rad';
end
if size(RA1,2)>1
    RA1  = celestial.coo.convertdms(RA1,'gH','r');
end
if size(Dec1,2)>1
    Dec1 = celestial.coo.convertdms(Dec1,'gD','R');
end
if size(RA2,2)>1
    RA2  = celestial.coo.convertdms(RA2,'gH','r');
end
if size(Dec2,2)>1
    Dec2 = celestial.coo.convertdms(Dec2,'gD','R');
end

switch lower(Units)
 case 'rad'
    % do nothing
 case 'deg'
    RA1   = RA1./RAD;
    Dec1  = Dec1./RAD;
    RA2   = RA2./RAD;
    Dec2  = Dec2./RAD;
 otherwise
    error('Unknown Units option');
end


% calculate Dist and PA:
[Dist,PA] = celestial.coo.sphere_dist(RA1,Dec1,RA2,Dec2);
if (Dist==0)
   OffsetLong = 0;
   OffsetLat  = 0;
else

   switch lower(OffsetType)
    case 'rd'
       % move from 1st to 2nd point first
       % in longitude then in latitude
   
       CD1  = celestial.coo.cosined([RA1,Dec1]);
       CD12 = celestial.coo.cosined([RA2,Dec1]);
       CD2  = celestial.coo.cosined([RA2,Dec2]);
   
       OffsetLong = acos(dot(CD1 ,CD12,2));
       OffsetLat  = acos(dot(CD12,CD2 ,2));
   
    case 'dr'
       % move from 1st to 2nd point first
       % in latitude then in longitude
   
       CD1  = celestial.coo.cosined([RA1,Dec1]);
       CD12 = celestial.coo.cosined([RA1,Dec2]);
       CD2  = celestial.coo.cosined([RA2,Dec2]);
   
       OffsetLat  = acos(dot(CD1 ,CD12,2));
       OffsetLong = acos(dot(CD12,CD2 ,2));
   
    otherwise
       error('Unknown OffsetType option');
   end
   
   OffsetLong = abs(OffsetLong).*sign(sin(PA));
   OffsetLat  = abs(OffsetLat) .*sign(cos(PA));

end

%OffsetLong.*RAD.*3600
%OffsetLat.*RAD.*3600
