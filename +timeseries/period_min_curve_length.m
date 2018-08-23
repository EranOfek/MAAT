function [pl,m]=period_min_curve_length(x,mn,mx,int,bin,c_x,c_y)
% Periodicity search using minimum curve length
% Package: timeseries
% Description: Search for periodicity in a time series, using the
%              minimum-curve-length method. The program calculates the
%              curve length for each trail frequency, and return the curve
%              length as function of frequency.
% Input  : - Data matrix, sorted by time.
%          - Minimum frequency to search.
%          - Maximum frequency to search.
%          - Frequency interval to search. defualt is 0.2/Time_Span.
%          - Optional binning. If 0, then don't use binning, if ~=0
%            then use bins of size 'bin'. Default is 0.
%          - The time column, defualt is 1.
%          - The dependent variable column, defualt is 2.
% Output : - Curve length as function of frequency,
%            [Frequency, Lengh].
%          - Frequency for which minimum length was obtained.
% Tested : Matlab 4.2
%     By : Eran O. Ofek                    Nov 1993
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

if nargin==3,
   c_x = 1;
   c_y = 2;
   bin = 0;
   int = 0.2./(max(x(:,c_x)) - min(x(:,c_x)))
elseif nargin==4,
   c_x = 1;
   c_y = 2;
   bin = 0;
elseif nargin==5,
   c_x = 1;
   c_y = 2;
elseif nargin==6,
   c_y = 2;
elseif nargin==7,
   % do nothing
else
   error('Illegal number of input arguments');
end


if bin>0.9,
   error('bin>0.9, bin is in phase');
end
j = 1;
for freq=mn:int:mx,
   f = folding(x,1./freq,c_x);
   if bin~=0,
      [f,tem1,tem2] = bining(f,bin,c_x);
   end
   pl(j,1) = k;
   pl(j,2) = curvlen(f);  %,c_x,c_y);
   k = k + int;
   j = j + 1;
end
[ml,ind] = min(pl(:,2:2));
m = pl(ind,1);
