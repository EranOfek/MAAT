function [X,N,Edges,Bounds]=realhist(V,Range,Option,Norm)
% Calculate histogram for a dataset in a given range.
% Package: Util.stat
% Description: Calculate histogram for a dataset in a given range.
% Input  : - Column vector for which to calculate the histogram.
%          - Histogram range, few optional formats:
%            [from, to, NumberOfBins] : for the histogram.
%                         (default is [min, max, 10]).
%            [X1;X2;...;Xn] : hist edges (see: histc_debug).
%          - Options :
%            'n'  : non cumulative (normal).
%            'c+' : cumulative histogram (in increasing direction).
%            'c-' : cumulative histogram (in decreasing direction).
%          - normalization -
%            'a' : as is (default).
%            'p' : probability normailzation. (normaliz sum to 1).
%            's' : spheric area normalization (range should be in radians).
%                  (the normalization is | sin(lat1) - sin(lat2) |).
%                  and the total sum is normalized to 1.
% Output : - Vector of histogram bin centeres.
%	   - Vector of number of data points in each bin.
%          - Edges vector, [X1;X2;...;Xn].
%          - One sigma Poisson's confidence interval per bin.
%            This is two column vector, for right and
%            left boundry confidence interval.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Sep 1995
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: R=rand(1000,1);
%          [x,n]=realhist(R,[0 1 50]);
%          bar(x,n);
% Reliable: 2
%--------------------------------------------------------------------
Eps     = 1e-8;
NBinDef = 10;
OneSig  = 0.6827;
if (nargin==1),
   MinVal = min(V(:,1));
   MaxVal = max(V(:,1));
   NormEps= (MaxVal-MinVal).*Eps;
   Range  = [MinVal-NormEps, MaxVal+NormEps, NBinDef];
   Option = 'n';
   Norm   = 'a';
elseif (nargin==2),
   Option = 'n';
   Norm   = 'a';
elseif (nargin==3),
   Norm   = 'a';
elseif (nargin==4),
   % no default
else
   error('Illegal number of input arguments');
end

if (length(Range(:,1))==1),
   % Range case
   if (length(Range)==2),
      Range(3) = NBinDef;
   end
   Step = (Range(2)-Range(1))./Range(3);
   X = Range(1) + Step.*0.5 + Step.*[0:1:Range(3)-1].';
   N = hist(V,X);
   N = N.';
   Edges = [Range(1):Step:Range(2)];
   EL = length(Edges);
else
   % Edges case
   % bin centers
   EL = length(Range);
   X  = 0.5.*(Range(1:EL-1) + Range(2:EL));
   N = histc_debug(V,Range);

   N = N(1:EL-1);
   Edges = Range;
end



% options
if (Option=='n'),
   % do nothing, default.
elseif (Option=='c+'),
   % cumulative (increasing)
   N = cumsum(N);
elseif (Option=='c-'),
   % cumulative (decreasing)
   CN = cumsum(rot90(N,2));
   N = rot90(CN,2);
else
   error('Unknown option');
end


if (nargout==4),
   % calculate Poisson errors
   Bounds = zeros(length(N),2);
   for I=1:1:length(N),
      [M, B] = poissfit(N(I),1-OneSig);
      Bounds(I,1:2) = B.';
   end
end



%normalizations
if (Norm=='a'),
   % default
elseif (Norm=='p'),
   SumN = sum(N);
   N = N./SumN;
   if (nargout==4),
      Bounds = Bounds./SumN;
   end
elseif (Norm=='s'),
   N = N./abs(sin(Edges(1:EL-1)) - sin(Edges(2:EL)));
   if (nargout==4),
      Bounds = Bounds./abs(sin(Edges(1:EL-1)) - sin(Edges(2:EL)));
   end
else
   error('Unknown normalization');
end


