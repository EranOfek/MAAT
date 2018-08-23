function Len=curvlen(Data)
% Calculate the length of a curve numerically.
% Package: Util.Geom
% Description: Calculate the length of a curve by summing the
%              distances (sqrt[X^2+Y^2]) between successive points.
% Input  : - Data matrix [X, Y], sorted by X.
% Output : - Curve length.
% Tested : Matlab 4.2
%     By : Eran O. Ofek                    Nov 1993
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Len=Util.Geom.curvlen(sortrows(rand(10,2),2));
% Reliable: 2
%--------------------------------------------------------------------------
ColX = 1;
ColY = 2;
%N    = size(Data,1);

DX = diff(Data(:,ColX));
DY = diff(Data(:,ColY));
Len = sum(DX.^2 + DY.^2);
