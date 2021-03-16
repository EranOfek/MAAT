function [RA1,RA2,Dec1,Dec2]=coo2box(RA,Dec,HalfSize,OutUnits)
% Calculate box vertices around coordinates (OBSOLETE: use coo2box)
% Package: celestial
% Description: Given a list of RA/Dec coordinates, and box half size,
%              calculate the approximnate positions of the box vertices
%              around the coordinates. Do not fix RA/Dec jumps.
% Input  : - Set of RA [radians].
%          - Set of Dec [radians].
%          - Box half size [RA, Dec] in radians.
%          - Output units {'rad','deg'}. Default is 'rad'.
% Output : - RA1 [radians]
%          - RA2 [radians]
%          - Dec1 [radians]
%          - Dec2 [radians]
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.coo.coo2box(1,1,0.01)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==3)
    OutUnits = 'rad';
end

HalfSize = HalfSize.*ones(1,2);

Dec  = Dec(:);
RA   = RA(:);

Dec1   = Dec - HalfSize(2);
Dec2   = Dec + HalfSize(2);
MaxDec = max(abs(Dec1),abs(Dec2));
RA1    = RA  - HalfSize(1)./cos(MaxDec);
RA2    = RA  + HalfSize(1)./cos(MaxDec);

RA1  = convert.angular('rad',OutUnits,RA1);
RA2  = convert.angular('rad',OutUnits,RA2);
Dec1 = convert.angular('rad',OutUnits,Dec1);
Dec2 = convert.angular('rad',OutUnits,Dec2);