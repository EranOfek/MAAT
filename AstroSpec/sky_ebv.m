function [Ebv,A]=sky_ebv(RA,Dec,CooType,Rv,CorrectHigh)
%--------------------------------------------------------------------------
% sky_ebv function                                               AstroSpec
% Description: Calculate the extinction from a local copy of the
%              Schlegel, Finkbeiner & Davis (1998) extinction maps.
% Input  : - J2000.0 RA or longitude [H M S] or [RAD] or sexagesimal.
%          - J2000.0 Dec or latitude [Sign D M S] or [RAD] or sexagesimal.
%          - Coordinates type:
%            'eq' : J2000.0 equatorial (default).
%            'g'  : Galactic.
%            'ec'  : ecliptic.
%          - R_{V}, default is 3.08 (for calculating A).
%          - {'y' | 'n'} Correct Schlegel et al. E(B-V) when >0.1,
%            using the Adams et al. (2013) correction.
%            I.e., at large E(B-V) values Schlegel et al. is probably
%            overestimating the extinction.
%            Default is 'y'.
% Output : - E(B-V) [mag].
%            The E(B-V) is returned from the point in the map which is
%            the nearest to the input coordinates.
%          - [A_U, A_B, A_V, A_R, A_I, A_J, A_H, A_K] (assuming R=3.08).
% Reference: Schlegel, Finkbeiner & Davis (1998; ApJ 500, 525).
% Needed : coco.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : hhtp://weizmann.ac.il/home/eofek/matlab/
% Example: [Ebv,A]=sky_ebv(1,1);
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;
Rv  = 3.08;
MapSize = 4096;
StepSize = 10000;

if (nargin==2),
   CooType = 'eq';
   Rv      = 3.08;
   CorrectHigh = 'y';
elseif (nargin==3),
   Rv      = 3.08;
   CorrectHigh = 'y';
elseif (nargin==4),
   CorrectHigh = 'y';
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end


RA  = convertdms(RA,'gH','r');
Dec = convertdms(Dec,'gD','r');

switch CooType
 case 'eq'
    % convert Equatorial to Galactic
    Coo = coco([RA, Dec],'j2000.0','g');
 case 'g'
    % do nothing - already in Galactic
    Coo = [RA, Dec];
 case 'ec'
    % convert Ecliptic to galactic
    Coo = coco([RA, Dec],'e','g');
 otherwise
    error('Unknown CooType Option');
end

L = Coo(:,1);    % Galactic longitude
B = Coo(:,2);    % Galactic latitude


%X = 2048.*sqrt(1 - sign(B).*sin(B)).*cos(L) + 2047.5;
%Y = -2048.*sign(B).*sqrt(1 - sign(B).*sin(B)).*sin(L) + 2047.5;
% correction to avoid values below 1 or above 4096...
X = 2047.5.*sqrt(1 - sign(B).*sin(B)).*cos(L) + 2048;
Y = -2047.5.*sign(B).*sqrt(1 - sign(B).*sin(B)).*sin(L) + 2048;
X = round(X);
Y = round(Y);

Ebv = zeros(size(L));

In = find(B<0);
Ip = find(B>=0);

if (isempty(In)==0),
   Im = fitsread('SFD_dust_4096_sgp.fits');
   for I=1:1:length(In),
      Ebv(In(I)) = Im(Y(In(I)),X(In(I)));
   end
end
if (isempty(Ip)==0),
   Im = fitsread('SFD_dust_4096_ngp.fits');
   for I=1:1:length(Ip),
      Ebv(Ip(I)) = Im(Y(Ip(I)),X(Ip(I)));
   end
end

% The Eb-v is probably overestimated when >0.1 mag
% Stanek 1998; Arce & Goodman 1999
% use correction from Adams et al. (2013):
switch lower(CorrectHigh)
 case 'y'
    Il = find(Ebv>0.1);
    Ebv(Il) = 0.1 + 0.65.*(Ebv(Il)-0.1);
 otherwise
    % do nothing
end

if (nargout>1),
   A_U = optical_extinction(Ebv,'B','V','U','C',Rv);
   A_B = optical_extinction(Ebv,'B','V','B','C',Rv);
   A_V = optical_extinction(Ebv,'B','V','V','C',Rv);
   A_R = optical_extinction(Ebv,'B','V','R','C',Rv);
   A_I = optical_extinction(Ebv,'B','V','I','C',Rv);
   A_J = optical_extinction(Ebv,'B','V','J','C',Rv);
   A_H = optical_extinction(Ebv,'B','V','H','C',Rv);
   A_K = optical_extinction(Ebv,'B','V','K','C',Rv);

   A = [A_U, A_B, A_V, A_R, A_I, A_J, A_H, A_K];
end



