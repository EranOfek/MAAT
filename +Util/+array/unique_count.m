function [UnVal,Count]=unique_count(Vec)
% Unique values and count the number of apperances of each value.
% Package: Util.array
% Description: Select unique values in numeric vector and count the
%              number of apperances of each value.
% Input  : - Numeric vector.
% Output : - Unique values
%          - Count of appearances per unique value.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [UnVal,Count]=Util.array.unique_count(Vec)
% Reliable: 2
%--------------------------------------------------------------------------

UnVal = unique(Vec);
Nun   = numel(UnVal);
Count = zeros(Nun,1);
for I=1:1:Nun
    Count(I) = sum(UnVal(I)==Vec);
end