function [R]=refraction(Alt,MetoData,Formula)
% Estimate atmospheric refraction, in visible light.
% Package: celestial.coo
% Description: Estimate atmospheric refraction, in visible light.
% Input  : - Vector of altitude [radians].
%          - [T (C), P (mb), Rh (%)]. default is [20, 1000, 0.8]
%            Relative humidity (Rh) is used only in 'Sa' formula.
%          - Formula type:
%            'AA' - Astronomical almanach, low accuracy for Alt>15deg,
%                   default.
%            'Sa' - Saastamoinen's formula for Alt>20deg.
% Output : - Refraction correction in radians. (Add to Alt).
% Reference : AA 2001
% Tested : Matlab 5.3
%     By : Eran O. Ofek                   Nov 2000
%    URL : http://weizmann.ac.il/home/eran/matlab/
% See also: refraction_wave.m
% Example: [R]=celestial.coo.refraction([20;30]./RAD);
% Reliable: 1
%------------------------------------------------------------------------------
RAD = 180./pi;

if (nargin==1)
   MetoData = [20, 1000, 0.8];
   Formula  = 'AA';
elseif (nargin==2)
   Formula  = 'AA';
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end
AltD = Alt.*RAD;

T = MetoData(1) + 273;  % convert to K
P = MetoData(2);

switch Formula
 case 'AA'
    Ia = find(AltD>15);
    Ib = find(AltD<=15 & AltD>-0.6);
    In = find(AltD<=-0.6);
    
    R = zeros(size(Alt));
   % R(Ia) = 0.00452.*P./(T.*tan(Alt(Ia)).*RAD);
    R(Ia) = 0.00452.*P./(T.*tand(AltD(Ia)+7.32./(AltD(Ia)+4.32)).*RAD);
    R(Ib) = P.*(0.1594 + 0.0196.*AltD(Ib) + 0.00002.*AltD(Ib).*AltD(Ib))./(T.*(1 + 0.505.*AltD(Ib) + 0.0845.*AltD(Ib).*AltD(Ib)).*RAD);
    R(In) = 0;

 case 'Sa'
    Rh    = MetoData(3);
    Delta = 18.36;
    Pw0   = Rh.*(T./247.1)^Delta;

    Q     = (P - 0.156.*Pw0)./T;

    Z     = pi./2 - Alt;
    TanZ  = tan(Z);
    R     = 16.271.*Q.*TanZ.*(1 + 0.0000394.*Q.*TanZ.^2) - 0.0000749.*P.*TanZ.*(1 + TanZ.^2);
    R     = R./3600./RAD;

 otherwise
    error('Unknown refraction formula option');
end
