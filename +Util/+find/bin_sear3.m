function IndVal=bin_sear3(X,Val)
%--------------------------------------------------------------------------
% bin_sear function                                                 FitFun
% Description: Binary search for a value in a sorted vector.
%              If the value does not exist, return the closes index.
% Input  : - sorted vector (ascending).
%	   - value to search.
% Output : - index of closest value.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    September 1994
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: IndVal=bin_sear([1:1:12]',5.4);
%--------------------------------------------------------------------------
N      = length(X);

CritN = 1e3;

if (N>CritN),
   VecInd       = (1:100:N).';
   [~,MinInd]   = min(abs(X(VecInd) - Val));
   VecInd       = (MinInd-100:MinInd+100).';
   VecInd       = VecInd(VecInd>0 & VecInd<=N);
   [~,MinInd]   = min(abs(X(VecInd) - Val));
   IndVal = VecInd(MinInd);
else
    VecInd = (1:1:N).';
    [~,IndVal] = min(abs(X(VecInd) - Val));
end