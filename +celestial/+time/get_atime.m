function T=get_atime(Date,GeodLong,UT1mUTC)
% Get current time, date, JD and LST.
% Package: celestial.time
% Description: Get current time, date, JD and LST.
% Input  : - Column vector of date in [JD], 
%            or [D M Y H M S] formats and in UTC time system.
%            If empty matrix then use current time and date.
%          - Geodetic east longitude [radians].
%          - UT1-UTC [s], default is 0.
% Output : - Structure of astronomical times, contains the following
%            fields:
%            .JD     - Julian day (UTC)
%            .Day    - Day in month (UTC)
%            .Month  - Month (UTC)
%            .Year   - Year (UTC)
%            .Hour   - Hour (UTC)
%            .Min    - Minutes (UTC)
%            .Sec    - Seconds (UTC)
%            .Frac   - Fraction of day (UTC)
%            .LST    - Local mean sidereal time [fraction of day]
%            .ISO    - String containing UTC date and time in
%                      standard ISO format.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jun 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: T=celestial.time.get_atime(convert.date2jd,35./RAD);
% Reliable: 1
%----------------------------------------------------------------------------- 

if (nargin<=2)
   UT1mUTC = 0;
end

if (isempty(Date)==1)
   Time = clock;
   T.Year   = Time(1);
   T.Month  = Time(2);
   T.Day    = Time(3);
   T.Hour   = Time(4);
   T.Min    = Time(5);
   T.Sec    = Time(6);
   T.JD     = convert.date2jd([T.Day, T.Month, T.Year, T.Hour, T.Min, T.Sec]);

%%% warning('get_atime: set time to -18.6 hours');
%%% T.JD = T.JD - 18.1./24;

else
   if (size(Date,2)==1)
      T.JD = Date;
   else
      T.JD = convert.date2jd(Date);
   end
   Time     = convert.jd2date(T.JD);

   T.Year   = Time(3);
   T.Month  = Time(2);
   T.Day    = Time(1);
   T.Frac   = Time(4);
   T_HMS    = convertdms(T.Frac,'f','H');
   T.Hour   = T_HMS(1);
   T.Min    = T_HMS(2);
   T.Sec    = T_HMS(3);
end

T.Frac = celestial.coo.convertdms([T.Hour, T.Min, T.Sec],'H','f');

T.LST  = celestial.time.lst(T.JD + UT1mUTC./86400,GeodLong,'m');

% Date and time in standard ISO format:
T.ISO  = sprintf('%04d-%02d-%02dT%02d:%02d:%06.3f',T.Year,T.Month,T.Day,T.Hour,T.Min,T.Sec);
