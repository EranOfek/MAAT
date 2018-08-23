function jd=julday1(Date)
% Convert Gregorian date in the range 1901 to 2099 to JD
% Package: celestial.time
% Description: Convert Gregorian date in the range 1901 to 2099 to
%              Julian days (see also: julday.m).
%              See also: convert.date2jd, celestial.time.julday
% Input  : - Gregorian date in the range 1901 to 2099 in one of the
%            following formats
%            [Day, Month, Year, Day_Fraction]
%            or [Day, Month, Year, Hour, Min, Sec]
%            or [Day, Month, Year] (in this case set Day_Fraction to 0.
% Output : Column vector of Julian Days.
% Tested : Matlab 3.5
%     By : Eran O. Ofek                   January 1994
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: celestial.time.julday1([1 1 2000 10 30 0]);
%          celestial.time.julday1([1 1 2000; 2 2 2000]);
% Reliable: 1
%--------------------------------------------------------------------------
y = Date(:,3);
m = Date(:,2);
d = Date(:,1);
switch size(Date,2)
 case 3
f = zeros(size(y));
 case 4
f = Date(:,4);
 case 6
 f = convertdms(Date(:,4:6),'H','f');
 otherwise
 error('Illegal number of column in Date matrix');
end

jd = 367.*y - floor(7.*(y+floor((m+9)./12))./4)+floor(275.*m./9)+d+1721013.5+f;
