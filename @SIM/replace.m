function Sim=replace(Sim,Range,Val,varargin)
% Replace pixels within specific values with another value.
% Package: @SIM
% Description: Given a SIM object image array, for each image replace any
%              value in the image with some specific range, or equal NaN,
%              with another value.
% Input  : - An array.
%          - List of ranges. This is a two columns matrix. Each row is
%            for a specific range, while the first column is for the lower
%            bound and the second column for the upper bound. If the first
%            column is NaN, then will search for NaN.
%          - A scalar. or a vector of the same length as the range matrix,
%            which values to use for the replacement.
% Output : - The array with the replaced elements.
% See also: Util.array.nan2val.m, replace.m SIM/nan2val.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim = replace(Sim,[-Inf 0;NaN 0],0);
% Reliable: 
%--------------------------------------------------------------------------

Sim = ufun2sim(Sim,@Util.array.replace,varargin{:},'FunAddPar',{Range,Val});

