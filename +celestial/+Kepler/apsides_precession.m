function DW=apsides_precession(Mass,A,E)
% First order estimation of the GR precession of the line of apsides
% Package: celestial.Kepler
% Description: First order estimation of the GR precession of the line
%              of apsides.
% Input  : - Central object mass [gram].
%          - Semi major axis [cm].
%          - Eccentricity.
% Output : - apsidel motion in radians per orbit.
% Reference: Danby 1989, Celestial Mechanics, p.68
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2002
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: DW=celestial.Kepler.apsides_precession(2e33,1e11,0.1)
% Reliable: 2
%------------------------------------------------------------------------------

SM = 1.9891e33;      % solar mass [gram]
AU = 1.49597870e13;  % au [cm]
G  = 6.672e-8;       % [cgs]
C  = 2.99792458e10;  % [cm]


% calculate the semi-latus rectum
P  = A.*(1-E.^2);

DW = 2.*pi.*3.*G.*Mass./(C.^2.*P);


% Barker & O'Connel (1975)

%(2.*pi./P).^(5./3) .*(G.*Msun./(c.^3)).^(2./3) .* Mc.*(4.*Mp + 3.*Mc)./(2.*(Mp + Mc).^(4./3)) ./(1-e.^2);
