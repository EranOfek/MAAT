function Nd=days_in_month(Year,Month)
% Return the number of days in month
% Package: celestial.time
% Description: Calculate the number of days in a given Gregorian or Julian
%              month.
% Input  : - Year (scalar).
%          - Month (scalar).
% Output : - Number of days in month. 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2003
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Nd=celestial.time.days_in_month(2000,2);
% Reliable: 2
%--------------------------------------------------------------------------

MonthE = Month + 1;
YearE  = Year;
if (MonthE==13),
   MonthE = 1;
   YearE  = Year + 1;
end

JD1 = convert.date2jd([1 Month  Year  0]);
JD2 = convert.date2jd([1 MonthE YearE 0]);

Nd = JD2 - JD1;
