function Mat=replace(Mat,Range,Val)
% Replace values in specific range, or equal NaN, with another value.
% Package: Util.array
% Description: Given an array replace any value in the array within some
%              specific range, or equal NaN, with another value.
% Input  : - An array.
%          - List of ranges. This is a two columns matrix. Each row is
%            for a specific range, while the first column is for the lower
%            bound and the second column for the upper bound. If the first
%            column is NaN, then will search for NaN.
%          - A scalar. or a vector of the same length as the range matrix,
%            which values to use for the replacement.
% Output : - The array with the replaced elements.
% See also: replace_with_noise.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Util.array.replace([1 2; NaN 2;17 19],[1 3],0);
%          Util.array.replace([1 2; NaN 2;17 19],[1 3;18 Inf],[0;NaN]);
% Reliable: 2
%--------------------------------------------------------------------------


Def.Val = 0;
if (nargin<2)
    Val = Def.Val;
end
if (isempty(Val))
    Val = Def.Val;
end

Nrange = size(Range,1);
Nval   = numel(Val);
for Irange=1:1:Nrange
    Ival = min(Irange,Nval);
    if (isnan(Range(Irange,1)))
        Mat(isnan(Mat)) = Val(Ival);
    else
        Mat(Mat>=Range(Irange,1) & Mat<=Range(Irange,2)) = Val(Ival);
    end
end

