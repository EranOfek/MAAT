function [PosConst]=constellation(Coo,Equinox)
% Find the constellations in which celestial coordinates are located.
% Package: celestial.stars
% Description: Find the constellations in which celestial coordinates are
%              located.
% Input  : - Matrix of [R.A., Dec.].
%            If two column matrix is given, then the first column is
%            the R.A. and the second is the Dec. (radians).
%            If Seven column matrix is given, then it is assumed
%            to be in sexagesimal format [H M S, Sign D M S].
%          - Equinox of coordinates, in JD (true eq. of date)
%            or 'J2000.0' for mean eq. of J2000.0.
%            'default is 'J2000.0'.
% Output : - Cell array with constellation names.
% Reference: Roman N.G., 1987, PASP 99, 695.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [PosConst]=celestial.stars.constellation([1 1;2 2]);
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;
ConstFileName = 'constellation.dat';

if (nargin==1),
   Equinox = 'J2000.0';
elseif (nargin==2),
   % do nothing
else
   error('Illigal number of input arguments');
end

SizeCoo = size(Coo);
N       = SizeCoo(1);
ColN    = SizeCoo(2);
if (ColN==2),
   % do nothing
elseif (ColN==7),
   RA  = celestial.coo.convertdms(Coo(:,1:3),'H','r');
   Dec = celestial.coo.convertdms(Coo(:,4:7),'D','R');
   Coo = [RA, Dec];
else
   error('Illigal number of columns in Coo matrix');
end

if (ischar(Equinox)==1),
   switch Equinox
    case 'J2000.0'
       Equinox = 2451545.0;
    otherwise
       error('Unkown equinox option');
   end
end

% precess coordinates to 1875.0
JD1875 = celestial.time.julday([1 1 1875 0]);

if (Equinox==2451545.0),
   % already in J2000.0
else
   % precess Coo to J2000.0
   Date   = celestial.time.jd2date(Equinox);
   JDY0   = celestial.time.julday([1 1 Date(3) 0]);
   Year   = Date(3) + (Equinox-JDY0)./365.25;
   EqYear = sprintf('j%06.1f',Year);
   Coo    = celestial.coo.coco(Coo,EqYear,'j2000.0');
end

% precess from mean equinox of J2000.0 to 1875.0
Coo = celestial.coo.coco(Coo,'j2000.0','j1875.0');
I   = find(Coo(:,1)<0);
Coo(I,1) = 2.*pi + Coo(I,1);

% read constellations file
%Const = cell(NLine,1);
%RA1   = zeros(NLine,1);
%RA2   = zeros(NLine,1);
%Dec1  = zeros(NLine,1);

FID = fopen(ConstFileName,'r');
%[RA1, RA2, Dec1, Const] = textscan(FID,'%f %f %f %s','CommentStyle','%');
[C] = textscan(FID,'%f %f %f %s','CommentStyle','%');
[RA1, RA2, Dec1, Const] = deal(C{1:4});
fclose(FID);

% convert to radians
RA1  = RA1.*15./RAD;
RA2  = RA2.*15./RAD;
Dec1 = Dec1./RAD;

% find constellation for each position
PosConst = cell(N,1);
for I=1:1:N,
   CI = find(Coo(I,1)>=RA1 & Coo(I,1)<RA2 & Coo(I,2)>Dec1);
   [~, MinInd] = min(Coo(I,2) - Dec1(CI));
   PosConst{I} = Const{CI(MinInd)};
   %Const{CI};  
end

