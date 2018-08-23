function [Par_E,Par_K]=kuiper_check(Date,RA,Dec,AbsK)
% Parallax due to Earth and object motion of a solar system object
% Package: celestial.SolarSys
% Description: Given a date and coordinates of a celestial object,
%              calculate the parallax due to the Earth motion (Par_E) 
%              and the parallax due to the object motion (Par_K), 
%              assuming the object is in e=0 orbit around the Sun,
%              its Sun distance is >>1au.
%              If (Par_E-Par_K)>Object_Proper_Motion then the object is
%              located outside the solar system.
% Input  : - Date [Day, Month, Year, fraction_of_day]
%            or [JD].
%          - RA [H M S] or [rad] or sexagesimal.
%          - Dec [Sign D M S] or [rad] or sexagesimal.
%          - sun-object distance [au], (default is 40au).
% Output : - Parallax due to Earth motion [arcsec/day].
%          - Parallax due to object motion [arcsec/day].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jun 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Par_E,Par_K]=celestial.SolarSys.kuiper_check([11 3 1999 0],[ 16 36 02],[+1 66 12 34]);
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;

if (nargin==3),
   AbsK = 40;  %au
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (size(Date,2)==4),
   JD = celestial.time.julday(Date);
else
   JD = Date;
end

RA  = celestial.coo.convertdms(RA,'gH','r');
Dec = celestial.coo.convertdms(Dec,'gD','R');

K = celestial.coo.cosined([RA,Dec]);
K = K.*AbsK;
K = K.';

[Coo,Vel]=celestial.SolarSys.calc_vsop87(JD,'Earth','e','E');


KEd  = K.'*Vel;
AbsV = sqrt(Vel.'*Vel);


CosPhi = KEd./(AbsK.*AbsV);
SinPhi = sqrt(1-CosPhi.^2);

%asin(SinPhi)*RAD;


Par_E = AbsV.*SinPhi./AbsK;

Par_K = 2.*pi./(365.25.*AbsK.^1.5);

Par_E = Par_E.*RAD.*3600;
Par_K = Par_K.*RAD.*3600;

