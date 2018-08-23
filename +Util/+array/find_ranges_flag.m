function [Flag]=find_ranges_flag(Vector,Ranges)
% Check if values are within some ranges.
% Package: Util.array
% Description: Given a vector and several ranges, return the a
%              vector indicating if a given position in the input vector
%              is included in one of the ranges.
% Input  : - Column vector of values.
%          - Matrix of ranges, in which each row specify a range.
%            The first column is for the range lower bound and the
%            second colum for the upper bound.
%            If the first column is NaN than look for NaN
% Output : - A vector (the same length as the input vector) in which
%            zero indicates that the position is not included in one
%            of the ranges, and other numbers indicate the index of
%            the range.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Flag]=find_ranges_flag((1:1:10),[1 3;6 7]);
%          [Flag]=find_ranges_flag([(1:1:10),NaN],[1 3;6 7; NaN NaN])
% Reliable: 1
%--------------------------------------------------------------------------

Nr = size(Ranges,1);
Flag = zeros(size(Vector));
for Ir=1:1:Nr
    if (isnan(Ranges(Ir,1)))
        Flag(isnan(Vector)) = Ir;
    else
        Flag((Vector>=Ranges(Ir,1) & Vector<=Ranges(Ir,2))) = Ir;
    end
end

