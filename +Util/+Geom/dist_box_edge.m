function MinDist=dist_box_edge(X,Y,EdgeX,EdgeY)
% Distance of points froma rectangular box.
% Package: Util.Geom
% Description: Given a rectangular box and a scalar position, calculate
%              the distance to the nearest rectangular edge.
% Input  : - Vector of X position.
%          - Vector of Y position.
%          - X position of edge [low, high].
%          - Y position of edge [low, high].
% Output : - Minimum distance.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: MinDist=Util.Geom.dist_box_edge([5;16],[15;3],[1 2048],[1 4096]);
% Reliable: 2
%--------------------------------------------------------------------------

MinX = min(abs(X-EdgeX(1)),abs(X-EdgeX(2)));
MinY = min(abs(Y-EdgeY(1)),abs(Y-EdgeY(2)));

MinDist = min(MinX,MinY);

