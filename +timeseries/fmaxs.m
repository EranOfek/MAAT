function [Res,Ind,WMeanPos]=fmaxs(Mat,ColX,ColY,MeanWindowSize);
%------------------------------------------------------------------------------
% fmaxs function                                                    timeseries
% Description: Given a matrix, find local maxima (in one of the columns) 
%              and return the maxima position and height.
% Input  : - Matrix of at least two columns.
%          - The column index of the independent variable,
%            default is 1.
%          - The column index of the dependent variable,
%            default is 2.
%          - Half window size for calculating the peak position
%            using wmean, default is 3.
% Output : - Two column matrix of all local maxima.
%            The first column is the independent variable
%            (maxima position) while the second column is the
%            dependent variable (maxima height).
%          - Vector of indices of the local maxima in the original
%            matrix.
%          - Weighted mean position, error and maximum of fitted parabola
%            for each peak.
%            The weights are taken as 1/sqrt(Height).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                  December 1993
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------------

if (nargin==1),
   ColX           = 1;
   ColY           = 2;
   MeanWindowSize = 3;
elseif (nargin==2),
   ColY           = 2;
   MeanWindowSize = 3;
elseif (nargin==3),
   MeanWindowSize = 3;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end

DiffVec = diff(sign(diff([0;Mat(:,ColY);0])));
Ind     = find(DiffVec==-2);

Res     = Mat(Ind,[ColX, ColY]);
LenMat  = size(Mat,1);

if (nargout>2),
   N  = length(Ind);   % number of peaks
   WMeanPos = zeros(N,3);
   for I=1:1:N,
      if ((Ind(I)-MeanWindowSize)<1 | (Ind(I)+MeanWindowSize)>LenMat),
         % can't calculate wmean, set to NaN
         Mean        = NaN;
         Err         = NaN;
         MaxParabola = NaN;
      else
         % Weighted mean
         SubMat     = Mat([Ind(I)-MeanWindowSize:1:Ind(I)+MeanWindowSize],:);
         [Mean,Err] = wmean([SubMat(:,1),1./sqrt(SubMat(:,2))]);

         % fit parabola
         Par = polyfit(SubMat(:,1),SubMat(:,2),2);
         MaxParabola = -Par(2)./(2.*Par(1));
      end
      WMeanPos(I,:) = [Mean, Err, MaxParabola];
   end
end
