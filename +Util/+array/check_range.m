function [NewMin,NewMax,MinChanged,MaxChanged]=check_range(DimLength,Min,Max)
% Replace out of bound indices with bound indices.
% Package: Util.array
% Description: Given the size of an N-dimensional array, and minimum and
%              maximum indices per each dimension, check if the indices
%              are out of bounds. If the minimum or maximum indices are
%              out of bounds return a new minimum and maximum indices
%              that are in bound and nearest to the original indices.
% Input  : - Vector of array size in each dimension (e.g., this can be
%            the output of size.m).
%          - Vector of minimum indices to check (one per dimension).
%          - Vector of maximum indices to check (one per dimension).
% Output : - Vector of new minimum indices.
%          - Vector of new maximum indices.
%          - Vector of flags indicating of the minimum indices were
%            changed (1) or not (0).
%          - Vector of flags indicating of the maximum indices were
%            changed (1) or not (0).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [NMin,NMax,MinCh,MaxCh]=check_range(size(rand(5,3,2)),[1 0 -1],[6 2 2]);
% Reliable: 2
%--------------------------------------------------------------------------
Ndim    = length(DimLength);
MinChanged = zeros(Ndim,1);
MaxChanged = zeros(Ndim,1);
NewMin = Min;
NewMax = Max;

% min
Is1 = find(Min<1);
NewMin(Is1)  = 1;
MinChanged(Is1) = 1;

% max
Ise = find(Max>DimLength);
NewMax(Ise) = DimLength(Ise);
MaxChanged(Ise) = 1;


