function [Mat,C,H]=contour_percentile(X,Y,Mat,Levels,varargin)
% Contour plot in which the contours are percentiles of the matrix sum.
% Package: plot
% Description: Given a matrix, generate a contour plot in which the contours
%              represent the percentiles of the total matrix sum. I.e., the
%              region encompassed within the contour contains the given
%              percentile of the matrix total sum.
% Input  : - Vector of X
%          - Vector of Y
%          - Matrix containing the number of counts in each cell.
%          - Percentiles levels, default is [0.6827 0.9545 0.9973].
%          - Optional line type (e.g., 'k-'), default is to use
%            contour.m defaults.
% Output : - Matrix containing the cumulative probability in each cell.
%          - Contour matrix as described in contourc.m
%          - Handle to a contourgroup object. This handle can
%            be used as input to clabel.m
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Jun 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mat,C,H]=plot.contour_percentile((1:1:10),(1:1:10),randi(100,10))
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==3)
   Levels   = 1-normcdf(0,[1;2;3],1).*2;
else
   % do nothing
end

Nlev = length(Levels);

Mat = Mat./Util.stat.sumnd(Mat);

PL = zeros(1,Nlev);
for Ilev=1:1:Nlev
   %PL(Ilev) = binsear_f('summatlevel',Levels(Ilev),[0 1],1e-3,Mat);
   PL(Ilev) = Util.find.fun_binsearch(@summatlevel,Levels(Ilev),[0 1],1e-3,Mat);
end

if (numel(varargin)==0)
   [C,H] = contour(X,Y,Mat,PL);
else
   [C,H] = contour(X,Y,Mat,PL,varargin{:});
end
