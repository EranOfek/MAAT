function [UT1mUTC,EOP]=ut1_utc(JD,Read)
% Return UT1-UTC (DUT1)
% Package: celestial.time
% Description: Return UT1-UTC (also known as DUT1).
% Input  : - Vector of Julian days (valid only after 1 1 1961).
%          - 'get' - get the latest EOP data from the IERS website and
%                    update the local version.
%            'use' - use local version of the EOP data (default).
% Output : - UT1-UTC [seconds].
%          - EOP table (see wget_eop.m for details).
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [UT1mUTC,EOP]=celestial.time.ut1_utc(2451545);
% Reliable: 2
%--------------------------------------------------------------------------


Def.Read = 'use';
if (nargin==1),
   Read = Def.Read;
end
InterpMethod = 'linear';

EOP = celestial.time.wget_eop(Read);

UT1mUTC = interp1(EOP.Cat(:,EOP.Col.MJD)+2400000.5,EOP.Cat(:,EOP.Col.UT1_UTC),JD,InterpMethod);
