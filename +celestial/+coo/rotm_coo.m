function RotM=rotm_coo(Type,EquinoxJD)
% Rotation matrix for coordinate conversion
% Package: celestial.coo
% Description: Generate a rotation matrix for coordinate conversion
%              and precession.
% Input  : - Type of rotation matrix
%            'e'  - equatorial with given equinox
%                   to ecliptic with the same ecliptic and equinox .
%                   Default is J2000.0
%            'E'  - ecliptic with given ecliptic and equinox 
%                   to the equatorial with the same equinox.
%                   Default is J2000.0
%            'G'  - galactic to equatorial with mean equinox of J2000.0
%            'g'  - equatorial with mean equinox of J2000.0 to galactic.
%            'p'  - precession matrix from mean equinox
%                   of date to mean equinox of J2000.0.
%            'P'  - precession matrix from mean equinox
%                   J2000.0 to mean equinox of date. 
%            'pd' - precession matrix from true equinox
%                   of date to mean equinox of J2000.0.
%            'Pd' - precession matrix from mean equinox
%                   J2000.0 to true equinox of date. 
%            'gSG'- Equatorial J2000.0 to Super Galactic.
%            'SGg'- Super Galactic to Equatorial J2000.0.
%            'gCMB'-Egalactic to WMAP 3rd year, CMB dipole.
%            'CMBg'- CMB dipole to Galactic.
%          - Equinox of coordinates (in Julian Day),
%            used only in the case of 'p' | 'P' | 'pd' | 'Pd' | 'ed' | 'Ed'
%            In case of 'E' or 'q' if this parameter is
%            not given it is taken as 2451545.0 (=J2000.0)
% Output : - Rotation matrix
% Reference : Ex. Supp. to the Astronomical Almanac.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: R=celestial.coo.rotm_coo('E');
%          % convert ecliptic coordinate [1 1] to equatorial:
%          celestial.coo.cosined([R*celestial.coo.cosined([1 1]).'].');
% Reliable: 1
%------------------------------------------------------------------------------

RADIAN = 180./pi;
J2000  = 2451545.5;

switch Type
 case {'e'}
    if (nargin==1)
       Obl  = celestial.coo.obliquity(J2000);
    else
       Obl  = celestial.coo.obliquity(EquinoxJD);
    end
    RotM = [1 0 0; 0 cos(Obl) sin(Obl); 0 -sin(Obl) cos(Obl)];
    %RotM = [-0.0548755604 +0.4941094279 -0.8676661490; -0.8734370902 -0.4448296300 -0.1980763734; -0.4838350155 +0.7469822445 +0.4559837762];
 case {'E','Ed'}
    if (nargin==1)
       Obl  = celestial.coo.obliquity(J2000);
    else
       Obl  = celestial.coo.obliquity(EquinoxJD);
    end
    RotM = [1 0 0; 0 cos(Obl) sin(Obl); 0 -sin(Obl) cos(Obl)].';
    %RotM = [-0.0548755604 +0.4941094279 -0.8676661490; -0.8734370902 -0.4448296300 -0.1980763734; -0.4838350155 +0.7469822445 +0.4559837762]';   
 case {'g'}
    RotM = [-0.0548755604 +0.4941094279 -0.8676661490; -0.8734370902 -0.4448296300 -0.1980763734; -0.4838350155 +0.7469822445 +0.4559837762]';
 case {'G'}
    RotM = [-0.0548755604 +0.4941094279 -0.8676661490; -0.8734370902 -0.4448296300 -0.1980763734; -0.4838350155 +0.7469822445 +0.4559837762];   
 case {'p','P','pd','Pd'}
    %T = (EquinoxJD - 2451545.0)./36525.0;
 
    [ZetaA,ZA,ThetaA]=celestial.coo.precession(EquinoxJD);
    
%     ZetaA  = 0.6406161.*T + 0.0000839.*T.*T + 0.0000050.*T.*T.*T;
%     ZA     = 0.6406161.*T + 0.0003041.*T.*T + 0.0000051.*T.*T.*T;
%     ThetaA = 0.5567530.*T - 0.0001185.*T.*T - 0.0000116.*T.*T.*T;
%     ZetaA  = ZetaA./RADIAN;
%     ZA     = ZA./RADIAN;
%     ThetaA = ThetaA./RADIAN;

    RotM = zeros(3,3);
    RotM(1,1) = cos(ZetaA).*cos(ThetaA).*cos(ZA) - sin(ZetaA).*sin(ZA);
    RotM(2,1) = cos(ZetaA).*cos(ThetaA).*sin(ZA) + sin(ZetaA).*cos(ZA);
    RotM(3,1) = cos(ZetaA).*sin(ThetaA);
    RotM(1,2) =-sin(ZetaA).*cos(ThetaA).*cos(ZA) - cos(ZetaA).*sin(ZA);
    RotM(2,2) =-sin(ZetaA).*cos(ThetaA).*sin(ZA) + cos(ZetaA).*cos(ZA);
    RotM(3,2) =-sin(ZetaA).*sin(ThetaA);
    RotM(1,3) =-sin(ThetaA).*cos(ZA);
    RotM(2,3) =-sin(ThetaA).*sin(ZA);
    RotM(3,3) = cos(ThetaA);
   
    switch Type
     case {'p'}
        RotM = RotM';
     case {'P'}
        RotM = RotM;
     case {'pd'}
        % calculate nutation matrix
        [~, NutMat]=celestial.coo.nutation(EquinoxJD);
        RotM = [NutMat']*[RotM'];
     case {'Pd'}
        % calculate nutation matrix
        [~, NutMat]=celestial.coo.nutation(EquinoxJD);
        RotM = NutMat*RotM;
     otherwise
        error('Unknown rotation matrix type');
    end    

 case {'gSG'}

    % Supergalactic coordinates. Supergalactic equator is conceptually
    % defined by the plane of the local (Virgo-Hydra-Centaurus)
    % supercluster, and the origin of supergalactic longitude is at the
    % intersection of the supergalactic and galactic planes. 
    % See the definition in RC2 (de Vaucouleurs et al. 1976) and RC3
    % (de Vaucouleurs et al. 1991): North SG pole at l=47.37 deg, b=6.32 deg.
    % Node at l=137.37, sgl=0 (inclination 83.68 deg).
    % Note: The position of the SG node reported in RC2 (137.29 deg) 
    % differs from RC3 (137.37 deg; see page 12, the difference is
    % mentioned but not commented).

    RotM = rotm(-83.68./RADIAN,'x')*rotm(-137.37./RADIAN,'z');

 case {'SGg'}

    RotM = [rotm(-83.68./RADIAN,'x')*rotm(-137.37./RADIAN,'z')].';

 case {'gCMB'}
    % Data from Table 6 in: Hinshaw et al. (2006)
    % 3.358\pm0.017 mK
    % l=263.86\pm0.04 deg
    % b=48.24\pm0.10

    RotM = rotm(48.24./RADIAN,'y')*rotm(-263.86./RADIAN,'z');

 case {'CMBg'}
    RotM = rotm_coo('gCMB').';

 otherwise
    error('Unknown rotation matrix type');
end
