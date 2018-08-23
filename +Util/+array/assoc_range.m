function RangeInd=assoc_range(Vec,Edges,Behav)
% Index of points in each bin.
% Package: Util.array
% Description: Given a vector of data points and a vector of edges,
%              for each value in the vector, return the index of the
%              bin (defined by the edges) to which it belongs.
%              
% Input  : - Vector od data points.
%          - Column vector of sorted edges.
%          - Behaviour type:
%            0 - Return NaN if the value is below the first edge or above
%                the last edge (default).
%            1 - Return 1 if the value is below the first edge, 2 if
%                the value is above the first edge and below the second
%                edge, etc.
% Output : - Indices of bins (defined by edges association for each
%            data point.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Vec = [-1;0;1.5;17]; Edges = [0 1 2 3].'; RangeInd=assoc_range(Vec,Edges)
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==2),
    Behav = 0;
end

RangeInd = zeros(size(Vec)).*NaN;
if (Behav==1),
   Edges = [-Inf; Edges; Inf];
end
Ned = numel(Edges);
for Ied=1:1:Ned-1,
    %Vec>Edges(Ied) & Vec<Edges(Ied+1)
    RangeInd(Vec>Edges(Ied) & Vec<=Edges(Ied+1)) = Ied;
end
