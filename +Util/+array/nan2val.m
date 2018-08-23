function Mat=nan2val(Mat,Val)
% Replace NaNs in an array with a specific values
% Package: Util.array
% Description: Replace NaNs in an array with a specific values
% Input  : - An array.
%          - A value which will replace all the NaN. Default is 0.
%            If empty then use default.
% Output : - The array with NaNs replaced.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Mat=nan2va([1 NaN;NaN 2],0);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Val = 0;
if (nargin<2)
    Val = Def.Val;
end
if (isempty(Val))
    Val = Def.Val;
end

Mat(isnan(Mat)) = Val;

