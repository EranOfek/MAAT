function Regions=flag2regions(Flags)
% Identify continous ranges of true values.
% Package: Util.array
% Description: Given a column vector flags (true/false) returns pairs of
%              indices of the positions of continuus regions in which the
%              flag is true.
% Input  : - Vector of flags (true/false).
% Output : - Two column matrix in which the eachj line indicates the
%            start and end index of a continuus region of true in the
%            input Flags vector.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Regions=flag2regions([0;1;1;1;0;1;1;1;0]);
% Reliable: 2
%--------------------------------------------------------------------------

Diff = diff([false; Flags; false]);
% region start with 1 and end with the next -1
RegionStart = find(Diff==1);
RegionEnd   = find(Diff==-1);
Regions     = [RegionStart, RegionEnd-1];

