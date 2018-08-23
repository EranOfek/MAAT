function [Ind]=find_ranges(Vector,Ranges)
% Indices of vector values found within one of several ranges.
% Package: Util.array
% Description: Given a vector and several ranges, return the indices of
%              values in the vector which are found within one of the
%              ranges.
% Input  : - Column vector of values.
%          - Matrix of ranges, in which each row specify a range.
%            The first column is for the range lower bound and the
%            second colum for the upper bound.
%            If the first column is NaN than look for NaN
% Output : - List of indices.
% Tested : Matlab 7.3
%     By : Eran O. Ofek               Feb 2008
%    URL : http://weizmann.ac.il/eofek/matlab/
% Example: [Ind]=Util.array.find_ranges([1 2 3 4 5 NaN],[0 1.5; NaN NaN]);
% Reliable: 1
%--------------------------------------------------------------------------

Nr = size(Ranges,1);
Ind = [];
for Ir=1:1:Nr
    if (isnan(Ranges(Ir,1)))
        Ind = [Ind; find(isnan(Vector))];
    else
        Ind = [Ind; find(Vector>=Ranges(Ir,1) & Vector<=Ranges(Ir,2))];
    end
end

