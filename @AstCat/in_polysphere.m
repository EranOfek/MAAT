function Flag=in_polysphere(AstC,Cols,Corners,Crit)
% Check if AstCat coordinates are within a spherical polygon.
% Package: @AstCat
% Description: Given an AstCat object with a single element check if
%              spherical coordinates in catalog are within a spherical
%              polygon.
% Input  : - AstCat object with a single element.
%          - Cell array of column names or vector of column indices
%            containing the Long/Lat info. Default is [1 2].
%            Default units of coordinates are radians.
%          - The verteces of a convex polygon on the celestial sphere
%            in whichits sides are great circles.
%            Each row correspond to one vertex.
%            This is a matrix of either 2 or 3 columns.
%            If two columns areprovided then these are [Long, Lat]
%            in radians. If three columns are given, these are
%            cosine directions.
%            The coordinates should be ordered such that the
%            right-hand rule is pointing toward the
%            center of the polygon.
%          - Flag indicating if to use ">" or ">=" in in_halfspace.m.
%            If 1 (default) then use (N dot R > C),
%            If 2 then use (N dot R) >= C
% Output : - A flag indicating if each row is in spherical polygon.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Flag=in_polysphere(AstC,[],Corners);
% Reliable: 
%--------------------------------------------------------------------------

if (isempty(Cols)),
    Cols = [1 2];
end
if (nargin<4),
    Crit = 1;
end

if (numel(AstC)>1),
    error('AstCat must contain a single element');
end

ColInd = colname2ind(AstC,Cols);
if (~isempty(AstC.ColUnits)),
    Units = AstC.ColUnits(ColInd(1));
    %ConvFactor = convert_units(Units,'rad');
    ConvFactor = convert.angular(Units,'rad');
else
    ConvFactor = 1;
end

Flag = in_polysphere(AstC.Cat(:,ColInd).*ConvFactor,Corners,Crit);


