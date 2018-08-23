function TI=thiele_innes(a,omega,Omega,I)
% Calculate the Thiele-Innes orbital elements
% Package: celestial.Kepler
% Description: Calculate the Thiele-Innes orbital elements.
% Input  : - a (semi major axis).
%          - omega (rad)
%          - Omega (rad)
%          - i (rad)
% Output : - Structure containing the 6 Thiele-Innes elements
%            (fields are A, B, C, F, G, H).
% See also: thiele_innes2el.m
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: TI=celestial.Kepler.thiele_innes(1,1,1,1);
% Reliable: 2
%--------------------------------------------------------------------------

TI.A = a.*( cos(omega).*cos(Omega) - sin(omega).*cos(I).*sin(Omega) );
TI.F = a.*(-sin(omega).*cos(Omega) - cos(omega).*cos(I).*sin(Omega) );
TI.B = a.*( cos(omega).*sin(Omega) + sin(omega).*cos(I).*cos(Omega) );
TI.G = a.*(-sin(omega).*sin(Omega) + cos(omega).*cos(I).*cos(Omega) );
TI.C = a.*sin(omega).*sin(I);
TI.H = a.*cos(omega).*sin(I);


% A + G = a.*(1+cos(i)).*cos(omega+Omega)
% A - G = a.*(1-cos(i)).*cos(omega-Omega)
% B - F = a.*(1+cos(i)).*sin(omega+Omega)
% -B- F = a.*(1-cos(i)).*sin(omega-Omega)

% (A+G).^2 + (B-F).^2 = a.^2.*(1+cos(i)).^2
% (A-G).^2 +(-B-F).^2 = a.^2.*(1-cos(i)).^2

% a = sqrt((A+G).^2 + (B-F).^2) + sqrt((A-G).^2 +(-B-F).^2)
% cos(i) = sqrt((A+G).^2 + (B-F).^2)./a - 1   % (two solutions)
% omega+Omega = atan2(A+G, B-F)
% omega-Omega = atan2(A-G, -B-F)




