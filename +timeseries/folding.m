function y=folding(x,p,c)
% Folding a timeseries by some period
  % Package: timeseries
% Description: Folding a time series into a period, and calculate the phase
%              for each data point in the time series.
% Input  : - Matrix in which one of the columns is the time.
%          - period to fold into.
%          - Column number to fold by (time column). Defualt is 1.
% Output : - Matrix similar to the input matrix, but in which
%            the time column is replaced by the phase., and
%            the time column is concatenate in the last column.
%            The output matrix is sorted by phase.
% Tested : Matlab 3.5
%     By : Eran O. Ofek                    Nov 1993
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: timeseries.folding(rand(100,2).*100,17);
% Reliable: 1
%--------------------------------------------------------------------------
if nargin==2,
   c=1;
elseif nargin==3,
   % do nothing
else
   error('Illegal number of input arguments');
end

jd_col      = length(x(1,:)) + 1;
y           = zeros(length(x(:,1)),jd_col);
TEMP        = x(:,c);
r           = TEMP;
TEMP        = TEMP./p-floor(TEMP./p);
y           = x;
y(:,c:c)    = TEMP;
y(:,jd_col) = r;
y           = sortrows(y,c);
