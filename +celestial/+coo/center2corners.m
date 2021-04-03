function [OutRA,OutDec]=center2corners(RA,Dec,H_FOV_RA,H_FOV_Dec)
% Return field corners given its center and size
% Package: celestial.coo
% Description: Given a field/s center and size calculate the field/s
%              corners.
% Input  : - Field center RA [radians].
%            This should be a column vector.
%          - Field center Dec [radians].
%            This should be a column vector.
%          - Half the field of view in RA [radians].
%          - Half the field of view in Dec [radians].
% Output : - Corners RA [radians] - 4 corners pear line.
%          - Corners Dec [radians] - 4 corners per line.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Mar 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [OutRA,OutDec]=celestial.coo.center2corners(1,1,0.01,0.01)
% Reliable: 2
%--------------------------------------------------------------------------

[OutRA1,OutDec1]=celestial.coo.add_offset(RA,Dec+H_FOV_Dec,H_FOV_RA,pi./2);
[OutRA2,OutDec2]=celestial.coo.add_offset(RA,Dec-H_FOV_Dec,H_FOV_RA,pi./2);
[OutRA3,OutDec3]=celestial.coo.add_offset(RA,Dec+H_FOV_Dec,H_FOV_RA,3.*pi./2);
[OutRA4,OutDec4]=celestial.coo.add_offset(RA,Dec-H_FOV_Dec,H_FOV_RA,3.*pi./2);

OutRA  = [OutRA1, OutRA2, OutRA3, OutRA4];
OutDec = [OutDec1, OutDec2, OutDec3, OutDec4];
