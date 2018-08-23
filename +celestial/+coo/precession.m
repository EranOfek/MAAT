function [ZetaA,ZA,ThetaA]=precession(JD)
% Calculate the Earth precession parameters
% Package: celestial.coo
% Description: Calculate the Earth precssion parameters as a function of
%              JD.
% Input  : - JD
% Output : - ZetaA [radians].
%          - ZA [radians]
%          - ThetaA [radians]
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Nov 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [ZetaA,ZA,ThetaA]=celestial.coo.precession(2451545+[0:1:5]');
% Reliable: 1
%------------------------------------------------------------------------------

InvRAD = pi./180;
T = (JD - 2451545.0)./36525.0;
 
ZetaA  = 0.6406161.*T + 0.0000839.*T.*T + 0.0000050.*T.*T.*T;
ZA     = 0.6406161.*T + 0.0003041.*T.*T + 0.0000051.*T.*T.*T;
ThetaA = 0.5567530.*T - 0.0001185.*T.*T - 0.0000116.*T.*T.*T;
ZetaA  = ZetaA.*InvRAD;
ZA     = ZA.*InvRAD;
ThetaA = ThetaA.*InvRAD;
