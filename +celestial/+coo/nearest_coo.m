function [NearCoo,NearInd,DistAll]=nearest_coo(CooList,Coo,Type)
%------------------------------------------------------------------------------
% nearest_coo function                                                   ephem
% Description: Given a list of coordinates (with arbitrary number of
%              dimensions), search for the coordinate in list which is
%              the nearest to a given (single) coordinate.
% Input  : - List of coordinates: [X] or [X Y] or [X Y Z],...
%          - Coordinate to search [X] or [X Y],...
%          - Coordinate type: {'plane' | 'sphere'}, default is 'plane' -
%            'sphere' is available only for 2-D [long, lat].
%            If 'sphere' then coordinates must be given in radians.
% Output : - Coordinate of nearest point in the list.
%          - Index of nearest point in the list.
%          - Distanace to all the points in list (in radians for 'sphere').
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Dec 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [NearCoo,NearInd,DistAll]=nearest_coo([1 1;2 2;3 1],[1 1]);
%          [NearCoo,NearInd,DistAll]=nearest_coo([1 1;2 2;3 1],[1 1.1],'sphere');
% Reliable: 2
%------------------------------------------------------------------------------
if (nargin==2),
   Type = 'plane';
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end


N = size(CooList,1);
switch Type
 case 'plane'
    %DistAll = sqrt(sum((CooList - ones(N,1)*Coo).^2,2));
    DistAll = sqrt(sum(bsxfun(@minus,CooList,Coo).^2,2));

 case 'sphere'
    DistAll = sphere_dist(Coo(1),Coo(2),CooList(:,1),CooList(:,2));

 otherwise
    error('Unknown Type Option');
end

[MinDist,NearInd] = min(DistAll);
NearCoo = CooList(NearInd,:);

