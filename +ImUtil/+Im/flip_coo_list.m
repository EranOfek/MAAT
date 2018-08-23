function Cat=flip_coo_list(Cat,Flip,Center,ColInd)
% Flip list of coordinates
% Package: ImUtil.Im
% Description: Flip the X and Y coordinates in a catalog around
%              a position. E.g., X=X-1024 or X=1024-X.
% Input  : - A catalog with [X,Y[ coordinates.
%          - A vector of two elements.
%            Elements are either -1 or 1. First element indicate if to
%            flip the X coordinate, and second is for the Y coordinate.
%            (e.g., If +1 then X=X-1024, if -1 then X=1024-X).
%          - A vector of two elements. 
%            Elements are the X and Y centers around the flip the
%            coordinates.
%          - Indices of X and Y columns.
%            Default is [1 2].
% Output : - AstCat object in which the X and Y columns are flipped.
% See also: AstCat/flip_center
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstC]=ImUtil.Im.flipp_coo_list(AstC,[-1 -1],[512 512],[1 2]);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<4)
    ColInd = [1 2];
end

if (Flip(1)<0)
    Cat(:,ColInd(1)) = Center(1) - Cat(:,ColInd(1));
else
    Cat(:,ColInd(1)) = Cat(:,ColInd(1)) - Center(1);
end

if (Flip(2)<0)
    Cat(:,ColInd(2)) = Center(2) - Cat(:,ColInd(2));
else
    Cat(:,ColInd(2)) = Cat(:,ColInd(2)) - Center(2);
end
    