function [OutRA,OutDec]=add_offset(RA,Dec,Offset,PA)
% Offset a position by angular distance and position angle
% Package: celestial.coo
% Description: Add an offset specified by angular distance and position
%              angle to celestial coordinates.
% Input  : - RA [radians].
%          - Dec [radians].
%          - Offset (angular distance) [radians].
%          - Position angle [radians].
% Output : - RA [radians].
%          - Dec [radians].
% See also: sphere_offset.m; sphere_dist.m
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Mar 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [OutRA,OutDec]=celestial.coo.add_offset(1,1,0.01,1)
% Reliable: 1
%--------------------------------------------------------------------------

[OutDec,OutRA]=reckon(Dec,RA,Offset,PA,'radians');
