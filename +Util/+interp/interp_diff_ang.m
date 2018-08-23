function [Yi]=interp_diff_ang(X,Y,Xi,Deg,Check)
% Stirling 4th order interpolation for angular values
% Package: Util.interp
% Description: Given a vector of time and a vector of coordinate on a sphere
%              interpolate the spherical coordinate in a list of times.
%              This function is talking into account
%              [0..2pi] angels discontinuity.
% Input  : - Equally spaced and (asendingly) sorted X.
%          - Y [angle in radians in range 0..2pi].
%          - X values for which to interpolate.
%          - degree of differences, default is 4.
%          - Check if X is equally spaced {'y' | 'n'}, default is 'n'.
% Output : - Interpolated Y values.
%            Return NaN in case that extrapolation is needed.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Seidelmann 1992, Explanatory Supp. to the Astron. Almanac
% See also: find_local_extramum.m; interp_diff.m
% Reliable: 2
%-----------------------------------------------------------------------------
DefDeg     = 4;
DefCheck   = 'n';
if (nargin==3)
   Deg    = DefDeg;
   Check  = DefCheck;
elseif (nargin==4)
   Check  = DefCheck;
elseif (nargin==5)
   % do nothing
else
   error('Illegal number of input arguments');
end


Y1 = sin(Y);
Y2 = cos(Y);

[Yi1] = Util.interp.interp_diff(X,Y1,Xi,Deg,Check);
[Yi2] = Util.interp.interp_diff(X,Y2,Xi,Deg,Check);

Yi    = atan2(Yi1,Yi2);

