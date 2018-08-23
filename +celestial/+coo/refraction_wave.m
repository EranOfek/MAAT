function [R,N]=refraction_wave(Alt,Lam,T,P,Pw)
%--------------------------------------------------------------------------
% refraction_wave function                                           ephem
% Description: Calculate the wavelength-dependent atmospheric refraction
%              and index of refraction based on Cox (1999) formula.
% Input  : - Altitude [radians].
%          - Wavelength [Ang]. Default is 5000 Ang.
%          - Temperature [C]. Default is 15 C .
%          - Pressure [hPa]. Default is 760 mm Hg.
%          - Partial vapour pressure. Default is 8 mm Hg.
% Output : - Atmospheric refraction [radians].
%          - Index of refraction.
% Reference: Filippenko 1982 (PASP 94, 715)
% Tested : Matlab R2011b
%     By : Eran O. Ofek                   Nov 2013
%    URL : http://weizmann.ac.il/home/eran/matlab/
% See also: refraction_wave.m
% Example: [R,N]=celestial.coo.refraction_wave([20;30]./RAD,7000);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Lam = 5000;
Def.T   = 15;
Def.P   = 760;
Def.Pw  = 8;

if (nargin==1)
    Lam   = Def.Lam;
    T     = Def.T;
    P     = Def.P;
    Pw    = Def.Pw;
elseif (nargin==2)
    T     = Def.T;
    P     = Def.P;
    Pw    = Def.Pw;
elseif (nargin==3)
        P     = Def.P;
    Pw    = Def.Pw;
elseif (nargin==4)
    Pw    = Def.Pw;
elseif (nargin==5)
    % do nothing
else
    error('Illegal number of input arguments');
end


Lam = Lam.*1e-4;  % Ang -> microns

% P = 760;
% T = 15;
% F = 8;

Nm1 = (64.328 + 29498.1./(146 - (1./Lam).^2) + 255.4./(41 - (1./Lam).^2))./1e6;
Nm1 = Nm1.*P.*(1 + (1.049 - 0.0157.*T).*1e-6.*P)./(720.883.*(1 + 0.003661.*T));
R1 = Nm1 - Pw.*(0.0624 - 0.000680./Lam.^2)./(1+0.003661.*T)./1e6;

N = R1 + 1;

% N = 1 + (6.4328e-5 + 2.94981e-2./(146-Lam.^-2) + 2.554e-4./(41-Lam.^-2)).* ...
%         (288.15./T).*(P./1013.25) - 4.349e-5.*(1-7.956.*1e-3.*Lam.^-2).* ...
%         (Pw./1013.25);
%     
%     
% R1 = N - 1; %(N.^2 - 2)./(2.*N.^2);

R  = R1.*tan(pi./2-Alt);

    