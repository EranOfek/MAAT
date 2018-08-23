function [TAImUTC,TTmUTC,TAI_UTC]=tai_utc(JD,Read)
% Get TAI-UTC (leap second)
% Package: celestial.time
% Description: Return the TAI-UTC time difference (leap second) for
%              a vector of Julian days. Also return TT-UTC.
% Input  : - Vector of JDs
%          - 'get' - get the latest TAI-UTC data from the IERS website and
%                    update the local version.
%            'use' - use local version of the TAI-UTC data (default).
%            See wget_tai_utc.m
% Output : - TAI-UTC time difference [seconds] for the requested JDs.
%            Return NaNs if out of range.
%          - TT-UTC [seconds].
%            Note that TT=TAI+32.184s.
%          - Data structure containing the TAI-UTC table.
%            See wget_tai_utc.m for additional information.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [TAImUTC,TTmUTC]=celestial.time.tai_utc([0;2451545;julday]);
% Reliable: 2
%--------------------------------------------------------------------------


Def.Read = 'use';
if (nargin==1)
   Read = Def.Read;
end

TAI_UTC = celestial.time.wget_tai_utc(Read);
Cat     = [[-Inf NaN NaN NaN]; TAI_UTC.Cat; [Inf NaN NaN NaN]];
Ind     = Util.array.assoc_range(JD,Cat(:,TAI_UTC.Col.JD));

TAImUTC = Cat(Ind,TAI_UTC.Col.DelTAU_UTC) + ...
          (JD - 2400000.5 - Cat(Ind,TAI_UTC.Col.T0)).*Cat(Ind,TAI_UTC.Col.S);

TTmUTC = 32.184 + TAImUTC;