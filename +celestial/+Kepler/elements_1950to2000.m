function [IN,OmegaN,omegaN]=elements_1950to2000(i,Omega,omega)
% B1950.0 FK4 orbital elements to J2000.0 FK5
% Package: celestial.Kepler
% Description: Convert solar system orbital elements given in the B1950.0
%              FK4 reference frame to the J2000.0 FK5 reference frame.
% Input  : - B1950.0 Orbital inclination [rad].
%          - B1950.0 Longitude of asecnding node [rad].
%          - B1950.0 Longitude of periastron [rad].
% Output : - J2000.0 Orbital inclination [rad].
%          - J2000.0 Longitude of asecnding node [rad].
%          - J2000.0 Longitude of periastron [rad].
% Reference: Seidelmann (1992), p. 314
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [IN,OmegaN,omegaN]=celestial.Kepler.elements_1950to2000(1,1,1)
% Reliable: 2
%--------------------------------------------------------------------------


RAD = 180./pi;
Lt = 4.50001688./RAD;
L  = 5.19856209./RAD;
J  = 0.00651966./RAD;

% B1950.0 to J2000.0
CosIN = cos(i).*cos(J) - sin(i).*sin(J).*cos(L+Omega);
IN      = acos(CosIN);
SinIN   = sqrt(1 - CosIN.^2);

omegaN  = omega + atan2(sin(J).*sin(L+Omega)./SinIN,...
                        (sin(i).*cos(J)+cos(i).*sin(J).*cos(L+Omega))./SinIN);

OmegaN  = -Lt + atan2(sin(i).*sin(L+Omega)./SinIN,...
                      (cos(i).*sin(J)+sin(i).*cos(J).*cos(L+Omega))./SinIN);
                  
                   
  
      




