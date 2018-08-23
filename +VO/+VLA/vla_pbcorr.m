function PBCorr=vla_pbcorr(Freq,Dist,InterpMethod)
% Calculate primary beam corrections for the VLA antena
% Package: VO.VLA
% Description: Calculate primary beam corrections for the VLA antena.
% Input  : - Frequency of observations [GHz].
%          - Angular distance from primary beam [arcmin].
%          - Interpolation method for frequency
%            {'linear'|'nearest'}, default is 'nearest'.
% Output : Primary beam correction factor (to multiply the flux by).
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jun 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

Def.InterpMethod = 'nearest';
if (nargin==2)
   InterpMethod = Def.InterpMethod;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end


Table = [  0.0738      -0.897  2.71   -0.242
           0.3275      -0.935  3.23   -0.378
           1.465       -1.343  6.579  -1.186
           4.885       -1.372  6.940  -1.309
           8.435       -1.306  6.253  -1.100
          14.965       -1.305  6.155  -1.030
          22.485       -1.417  7.332  -1.352
          43.315       -1.321  6.185  -0.983];
% distance [arcmin] above to ignore correction (set to NaN)
Cutoff = [599.6
          135.1
           29.8
            9.13
            5.25
            2.95
            1.97
	    1.02];

IFreq = interp1(Table(:,1),Table(:,1),Freq,InterpMethod);
PB3 = interp1(Table(:,1),Table(:,2),Freq,InterpMethod);
PB4 = interp1(Table(:,1),Table(:,3),Freq,InterpMethod);
PB5 = interp1(Table(:,1),Table(:,4),Freq,InterpMethod);
Cut = interp1(Table(:,1),Cutoff,Freq,InterpMethod);

DistF = (Dist.*IFreq).^2;

PBCorr = 1 + PB3.*1e-3.*DistF + ...
             PB4.*1e-7.*DistF.^2 + ...
             PB5.*1e-10.*DistF.^3;

PBCorr(Dist>Cut) = NaN;

