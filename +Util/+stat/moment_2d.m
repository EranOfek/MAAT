function [M1,M2]=moment_2d(Mat,X,Y)
%------------------------------------------------------------------------------
% moment_2d function                                                 AstroStat
% Description: Calculate first and second moments of a 2D matrix.
% Input  : - Matrix
%          - X coordinate of matrix j-axis. Default is [1:1:size(Mat,2)].'
%          - Y coordinate of matrix i-axis. Default is [1:1:size(Mat,1)].'
% Output : - First moment [E(x),E(y)].
%          - Second Moment [E(x^2) E(x*y); E(y*x) E(y^2)].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jun 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [M1,M2]=moment_2d(rand(10,10));
% Reliable: 2
%------------------------------------------------------------------------------

Def.X = [1:1:size(Mat,2)].';
Def.Y = [1:1:size(Mat,1)].';
if (nargin==1),
   X   = Def.X;
   Y   = Def.Y;
elseif (nargin==2),
   Y   = Def.Y;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end
   

[MatX,MatY] = meshgrid(X,Y);
SumMat      = sum(sum(Mat));

M1(1)   = sum(sum(Mat.*MatX))./SumMat;
M1(2)   = sum(sum(Mat.*MatY))./SumMat;

if (nargout>1),
   M2      = zeros(2,2);
   M2(1,1) = sum(sum(Mat.*(MatX-M1(1)).^2))./SumMat;
   M2(2,2) = sum(sum(Mat.*(MatY-M1(2)).^2))./SumMat;
   M2(1,2) = sum(sum(Mat.*(MatX-M1(1)).*(MatY-M1(2))))./SumMat;
   M2(2,1) = M2(1,2);
end
