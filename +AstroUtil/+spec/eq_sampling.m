function [NewList1,NewList2]=eq_sampling(List1,List2,Samp,Method)
% Resample two [X,Y] lists to have the same sampling (x).
% Package: AstroUtil.spec
% Description: Given two lists, each contains [X,Y], equalize the sampling
%              frequency of the two lists by interpolating both lists to
%              at a specified X.
% Input  : - First list, in which the first column is the independent
%            varaible, and the second column is the value for each point.
%            If this is a cell array, then the first cell contains
%            a column vector of the independent variable, and
%            the second cell contain a matrix of the dependt variable,
%            and the interpolation is preformed for each column separately.
%          - Second list, (like the first list).
%          - New sampling rate, if a scalar is then taken as a constant
%            sampling rate. In case that a vector is given, it is taken as
%            the new independent variable. In case that empty matrix ([])
%            is given then the minimum sampling rate of the two lists
%            is taken. (default).
%          - Interpolation method:
%            'linear'  - linear interpolation (default).
%            'nearest' - nearest neighbor interpolation
%            'spline'  - cubic spline interpolation
%            'cubic'   - cubic interpolation
% Output : - The first list with the new sampling. 
%            If the input was a cell arry then the output is a cell
%            array too. 
%          - The second list with the new sampling.
%            If the input was a cell arry then the output is a cell
%            array too. 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 2000   
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2)
   Samp   = [];
   Method = 'linear';
elseif (nargin==3)
   Method = 'linear';
elseif (nargin==4)
   % no default.
else
   error('Illegal number of input arguments');
end


if (iscell(List1)==1)
   List1X = List1{1};
   List1Y = List1{2};
else
   List1X = List1(:,1);
   List1Y = List1(:,2);
end

if (iscell(List2)==1)
   List2X = List2{1};
   List2Y = List2{2};
else
   List2X = List2(:,1);
   List2Y = List2(:,2);
end


MinVal = max([min(List1X), min(List2X)]);
MaxVal = min([max(List1X), max(List2X)]);
if (isempty(Samp))
   Rate = min([diff(List1X); diff(List2X)]);

   SampVec = [MinVal:Rate:MaxVal].';
else
   if (length(Samp)==1)
      SampVec = [MinVal:Samp:MaxVal].';
   else
      if strcmp(Method,'spline')||strcmp(Method,'cubic')
         % Those method support extrapolation
         SampVec = Samp;
      else
         % Truncate Samp to the lists common interval to avoid NaN result
         % outside this interval
         SampVec = Samp(Samp>=MinVal&Samp<=MaxVal);
      end
   end
end

Y1 = interp1(List1X, List1Y, SampVec, Method);
Y2 = interp1(List2X, List2Y, SampVec, Method);

if (iscell(List1)==1)
   NewList1{1} = SampVec;
   NewList1{2} = Y1;
else
   NewList1 = [SampVec(:), Y1(:)];
end

if (iscell(List2)==1)
   NewList2{1} = SampVec;
   NewList2{2} = Y2;
else
   NewList2 = [SampVec(:), Y2(:)];
end

