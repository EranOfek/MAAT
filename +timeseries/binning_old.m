function [ResMat,NewMat]=binning_old(Data,BinSize,Start,End)
% The old version of the binning function. Use binning instead.
% Package: timeseries
% Description: Binning a timeseries using equally spaced bins.
%              In each bin, calculate the mean, median, std,
%              skewness of the data.
%              OBSOLETE: Use binning.m instead.
% Input  : - Matrix, in which the first column is the
%            time, the second column is the observation, and
%            optional third column is the value error.
%          - bin size, in units of the "time" column.
%          - First bin start Time (default is start point).
%          - Last bin end Time (default is end point).
% Output : - Matrix with the following columns:
%            [Bin_Center,
%             <Y>,
%             StD(Y)/sqrt(N),
%             <X>,
%             Weighted-<Y>,
%             Formal-Error<Y>,
%             N,
%             median(Y),
%             StD(Y),
%             range(X)/2]
%            In case that no errors are given, columns 5 and 6
%            will be zeros.
%          - The matrix as above, but the NaNs are eliminated.
% See also: bin_by_eye.m, runmean.m, bin_en.m
% Tested : Matlab 5.0
%     By : Eran O. Ofek                    Feb 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%--------------------------------------------------------------------------
import Util.stat.*

Col.X = 1;
Col.Y = 2;
Col.E = 3;

Def.Start = min(Data(:,Col.X));
Def.End   = max(Data(:,Col.X));
if nargin==2,
   Start = Def.Start;
   End   = Def.End;
elseif nargin==3,
   End   = Def.End;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end

% if (length(Data(1,:))==2),
%    % without errors
%    Err = 'n';
% else
%    Err = 'y';
% end

%Ntot = length(Data(:,Col.X));

Nbin = floor((End-Start)./BinSize);

% Bin_Center; <X>; <Y>; <Y^2>; N
ResMat = zeros(Nbin,10).*NaN;

% initialize bin center
ResMat(:,1) = [(Start+0.5.*BinSize):BinSize:(End-0.5.*BinSize)]';

for Ib=1:1:Nbin,
   If = abs(ResMat(Ib,1)-Data(:,Col.X))<(0.5.*BinSize);
   Nf = sum(If);
   ResMat(Ib,2) = mean( Data(If,Col.Y) );
   ResMat(Ib,3) = std( Data(If,Col.Y) )./sqrt(Nf);
   ResMat(Ib,4) = mean( Data(If,Col.X) );
   if (size(Data,2)>=Col.E),
      [M,E] = wmean(Data(If,Col.Y),Data(If,Col.E));
      ResMat(Ib,5) = M;
      ResMat(Ib,6) = E;
   end
   ResMat(Ib,7) = Nf;
   ResMat(Ib,8) = median( Data(If,Col.Y) );
   ResMat(Ib,9) = std( Data(If,Col.Y) );
   if (Nf<2),
      ResMat(Ib,10) = 0;
   else
      ResMat(Ib,10) = range(Data(If,Col.X))./2;
   end
end


if (nargout>1),
   % eliminate NaNs
   Innan  = isnan(ResMat(:,2))==0;
   NewMat = ResMat(Innan,:);
end

   


