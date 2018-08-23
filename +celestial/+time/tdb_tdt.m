function Diff=tdb_tdt(JD)
% Approximate TDB-TT
% Package: celestial.time
% Description: Calculate approximate difference between TDT and TDB
%              time scales.
% Input  : - Vector of julian days.
% Output : - TDB-TDT (TDB-TT) [seconds].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jun 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Diff=celestial.time.tdb_tdt(julday([1 1 2000]));
% Reliable: 2
%-----------------------------------------------------------------------

RAD = 180./pi;
%if (nargin==1),
%   Option = 'b';
%end

G = 357.53 + 0.9856003.*(JD - 2451545.0);
G = G./RAD;
Diff = 0.001658.*sin(G) + 0.000014.*sin(2.*G);  %sec.

