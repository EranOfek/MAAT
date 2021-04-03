function AM=hardie(X,Algo)
% The Hardie airmass formula
% Package: celestial.coo
% Description: Calculate airmass using the Hardie formula.
% Input  : - Matrix of zenith distances [radians].
%          - Algorith: {'hardie','csc'}. Default is 'hardie'.
% Output : - Air mass.
% Tested : Matlab 4.2
%     By : Eran O. Ofek                    Jan 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: airmass.m; hardie_inv.m
% Example: AM=celestial.coo.hardie([1;1.1]);
% Reliable: 1
%--------------------------------------------------------------------------

RAD = 180./pi;

Def.Algo = 'hardie'; 
if nargin==1
   Algo = Def.Algo;
elseif nargin==2
   % do nothing
else
   error('Illegal number of input arguments');
end

SecZ = 1./cos(X);
switch lower(Algo)
    case 'hardie'
        AM   = SecZ - 0.0018167.*(SecZ - 1) - 0.002875.*(SecZ - 1).*(SecZ - 1) - 0.0008083.*(SecZ - 1).*(SecZ - 1).*(SecZ - 1);
    case 'csc'
        AM = SecZ;
    otherwise
        error('Unknown airmass algorith option');
end

FlagI = X>87./RAD;
AM(FlagI) = NaN;



