function [AstC]=flip_center(AstC,Flip,Center,ColNames)
%--------------------------------------------------------------------------
% flip_center function                                       class/@AstCat
% Description: Flip the X and Y coordinates in an AstCat object around
%              a position. E.g., X=X-1024 or X=1024-X.
% Input  : - AstCat object.
%          - Two column matrix. Row per element in the input AstCat object.
%            Elements are either -1 or 1. First column indicate if to
%            flip the X coordinate, and second is for the Y coordinate.
%            (e.g., If +1 then X=X-1024, if -1 then X=1024-X).
%          - Two column matrix. Row per element in the input AstCat object.
%            Elements are the X and Y centers around the flip the
%            coordinates.
%          - Cell array of X and Y column names or vector of column
%            indices. Default is {'XWIN_IMAGE','YWIN_IMAGE'}.
% Output : - AstCat object in which the X and Y columns are flipped.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstC]=flip_center(AstC,[-1 -1],[512 512],{'XWIN_IMAGE','YWIN_IMAGE'})
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<4)
    ColNames = {'XWIN_IMAGE','YWIN_IMAGE'};
end

Nc = numel(AstC);
Center = bsxfun(@times,Center,ones(Nc,1));
Flip   = bsxfun(@times,Flip,  ones(Nc,1));
for Ic=1:1:Nc
    % for each AstCat catalog
    
    ColInd = colname2ind(AstC(Ic),ColNames);
    
    if (Flip(Ic,1)<0)
        AstC(Ic).Cat(:,ColInd(1)) = Center(Ic,1) - AstC(Ic).Cat(:,ColInd(1));
    else
        AstC(Ic).Cat(:,ColInd(1)) = AstC(Ic).Cat(:,ColInd(1)) - Center(Ic,1);
    end
    
    if (Flip(Ic,2)<0)
        AstC(Ic).Cat(:,ColInd(2)) = Center(Ic,2) - AstC(Ic).Cat(:,ColInd(2));
    else
        AstC(Ic).Cat(:,ColInd(2)) = AstC(Ic).Cat(:,ColInd(2)) - Center(Ic,2);
    end
    
    AstC(Ic).SortedBy    = [];
    AstC(Ic).SortedByCol = [];
end
