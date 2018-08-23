function [RA3,Dec3]=pole_from2points(RA1,Dec1,RA2,Dec2)
% Find pole of a great circle defined by two points on the sphere.
% Package: celestial.coo
% Description: Given two points on the celestial sphere (in any system)
%              describing the equator of a coordinate system,
%              find one of the poles of this coordinate system.
% Input  : - RA of 1st point [radians].
%          - Dec of the 1st point [radians]
%          - RA of the 2nd point [radians]
%          - Dec of the 2nd point [radians]
% Output : - RA of one of the poles.
%          - Dec of the pole.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jul 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [RA3,Dec3]=celestial.coo.pole_from2points(RA1,Dec1,RA2,Dec2)
% Reliable: 2
%------------------------------------------------------------------------------


[~,PA] = celestial.coo.sphere_dist_fast(RA1,Dec1,RA2,Dec2);
Phi = (2.*pi-PA) - pi./2;
[Dec3,RA3] = reckon(Dec1,RA1,pi./2,Phi,'radians');



% CosD  = sin(Dec1).*sin(Dec2) + cos(Dec1).*cos(Dec2).*cos(RA2-RA1);
% D     = acos(CosD);
% Gamma = asin(cos(Dec2).*sin(RA2-RA1)./sin(D)) - pi./2;
% FlagD = Dec2<Dec1;
% Gamma(FlagD) = pi./2 - asin(cos(Dec2(FlagD)).*sin(RA2(FlagD)-RA1(FlagD))./sin(D(FlagD)));
% 
% 
% Dec3  = asin(cos(Dec1).*cos(Gamma));
% RA3   = RA1 + atan2(sin(Gamma)./cos(Dec3),-sin(Dec1).*cos(Gamma)./cos(Dec3));
