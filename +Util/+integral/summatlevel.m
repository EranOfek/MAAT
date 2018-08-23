function Frac=summatlevel(Level,Mat)
%--------------------------------------------------------------------------
% summatlevel function                                             General
% Description: Given a matrix and a level value, return the sum of all the
%              values in the matrix which are larger than the level.
% Input  : - Level.
%          - Matrix.
% Output : - Sum of values larger than level.
% See Also: contour_percentile.m
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Jun 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

I    = Mat>=Level;
Frac = sumnd(Mat(I));

