function [Ev,UBV]=stellar_tracks(Mass,Metal,Type)
% Geneva stellar tracks as a function of time.
% Package: AstroUtil.stars
% Description: Given an initial mass and metllicity return the Geneva
%              stellar tracks as a function of time.
% Input  : - Initial mass [solar mass].
%            Can be one of the followings:
%            [0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |
%             1.25 | 1.5 | 1.7 | 2 | 2.5 | 3 | 4 | 5 | 7 | 9 | 10 |
%             12 | 15 | 20 | 25 | 40 | 60 | 85 | 120].
%          - Metallicity [mass fraction].
%            Can be one of the followings:
%            [0.0004 | 0.001 | 0.004 | 0.008 | 0.02 | 0.04 | 0.1]
%          - Model type ["c" | "e" | "l" | "m" | "p" | "o"],
%            see table1.dat for details. Default is "c".
% Output : - Evolution structure, contains the following fields:
%            .Number - grid point number
%            .Age    - age [yr]
%            .Mass   - actual mass [solar mass]
%            .LogL   - log(Luminosity) [log(solar lum)]
%            .LogT   - log(Effective Temp.) [log(K)]
%            .H; .He; .C12; C13; .N14; .O16; .O17; .O18; .Ne20; .Ne22
%                    - surface abundance [mass fraction] of these elements.
%            .CoreT  - Stellar core temperature for Wolf-Rayet stars [K].
%            .ML     - Mass loss rate [solMass/yr].
%          - UBV structure, contains the following fields:
%            .Number - grid point number
%            .Age    - age [yr]
%            .Mass   - actual mass [solar mass]
%            .LogL   - log(Luminosity) [log(solar lum)]
%            .LogT   - log(Effective Temp.) [log(K)]
%            .Logg   - log(Surface gravity) [log(cm/s^2)]
%            .U; .B; .V; .R; .I; .J; .H; .K; .L; .Lt; .M
%                    - Abs. magnitude in the respective bands.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Feb 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Ev,UBV]=AstroUtil.stars.stellar_tracks(1,0.04,'c');
% Reliable: 2
%--------------------------------------------------------------------------
DefType = 'c';

if (nargin==2)
   Type  = DefType;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

% convert Type ('celmpo') to type index
TypeVec = 'celmpo';
TypeInd = strfind(TypeVec,Type);
if (isempty(TypeInd)==1)
   error('Unknown Type option');
end

%--- get Geneva model tracks ---
%-----------------
%--- evolution ---
%-----------------
%load evol.mat;
evol = cats.stars.GenevaTracks.evol;
I = find(evol(:,2)==Mass & evol(:,1)==Metal & evol(:,3)==TypeInd);
if (isempty(I)==1)
   error('Model parameters not available');
end
Ev.N      = evol(I,4);
Ev.Age    = evol(I,5);
Ev.Mass   = evol(I,6);
Ev.LogL   = evol(I,7);
Ev.LogT   = evol(I,8);
Ev.H      = evol(I,9);
Ev.He     = evol(I,10);
Ev.C12    = evol(I,11);
Ev.C13    = evol(I,12);
Ev.N14    = evol(I,13);
Ev.O16    = evol(I,14);
Ev.O17    = evol(I,15);
Ev.O18    = evol(I,16);
Ev.Ne20   = evol(I,17);
Ev.Ne22   = evol(I,18);
Ev.CoreT  = evol(I,19);
Ev.ML     = evol(I,20);


%-----------------
%--- evolution ---
%-----------------
%ubv=load2('ubv.mat');
ubv = cats.stars.GenevaTracks.ubv;
I = find(ubv(:,2)==Mass & ubv(:,1)==Metal & ubv(:,3)==TypeInd);
UBV.N      = ubv(I,4);
UBV.Age    = ubv(I,5);
UBV.Mass   = ubv(I,6);
UBV.LogL   = ubv(I,9);
UBV.LogT   = ubv(I,7);
UBV.Logg   = ubv(I,8);
UBV.V      = ubv(I,10);
UBV.B      = ubv(I,12) + UBV.V;
UBV.U      = ubv(I,11) - ubv(I,12) - UBV.V;
UBV.R      = UBV.V - ubv(I,13);
UBV.I      = UBV.V - ubv(I,14);
UBV.K      = UBV.V - ubv(I,15);
UBV.H      = UBV.K + ubv(I,19);
UBV.J      = UBV.H + ubv(I,18);
UBV.L      = UBV.K - ubv(I,20);
UBV.Lt     = UBV.J - ubv(I,23);
UBV.M      = UBV.K - ubv(I,24);
