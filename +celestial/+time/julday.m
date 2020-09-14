function JD=julday(Date,Output)
% Convert Julian/Gregorain date to JD
% Package: celestial.time
% Description: Convert Julian/Gregorian date to Julian Day.
%              OBSOLETE: Use convert.date2jd instead.
%              See also: convert.time
% Input  : - Gregorian of Julian date in one of the following formats
%            [Day, Month, Year, Day_Fraction]
%            or [Day, Month, Year, Hour, Min, Sec]
%            or [Day, Month, Year] - in this case set Day_Fraction to 0.
%            Alternatively this can be a string or a cell array of strings
%            in which each string contains the date in the format:
%            'yyyy-mm-ddTHH:MM:SS' (e.g., '2010-08-11T15:01:56') or:
%            'yyyy-mm-dd HH:MM:SS'.
%            If argument is not provided then the program will calculate
%            the JD for now using the clock UTC computer (time zone
%            included).
%          - Output type. Options are:
%            'JD'  - Julian days (default).
%            'MJD' - Modified JD (JD-2400000.5).
% Output : - Row vector of Julian days.
% Tested : Matlab 3.5
%     By : Eran O. Ofek                    Jan 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.time.julday([1 1 2000 10 30 0]);
%          celestial.time.julday([1 1 2000; 2 2 2000]);
%          celestial.time.julday('2010-10-12 10:10:10.111');
%          celestial.time.julday({'2010-10-12 10:10:10.111'});
%          celestial.time.julday;  % JD of now
% Reliable: 1
%--------------------------------------------------------------------------

if (nargin==0)
   Date = clock;
   Date = Date(:,[3 2 1 4 5 6]);
end
if (isempty(Date))
   %Date = clock;
   % read UTC time
   Date = datevec(datetime('now', 'TimeZone', 'UTC'));
   Date = Date(:,[3 2 1 4 5 6]);
end
if (ischar(Date) || iscell(Date))
    %Date=date_str2vec(Date);
    Date = convert.str2date (Date);
    SizeD = size(Date,2);
    if (SizeD==3)
        Date = Date(:,[3 2 1]);
    elseif (SizeD==6)
        Date = Date(:,[3 2 1 4 5 6]);
    else
        error('unknown number of columns in date format - support 3 or 6');
    end
end
%    Date = {Date};
% end
% if (iscell(Date)),
%     if (~strcmp(Date{1}(11),'T')),
%         Date = datevec(Date,'yyyy-mm-dd HH:MM:SS');
%     else
%         Date = datevec(Date,'yyyy-mm-ddTHH:MM:SS');
%     end
%     Date = Date(:,[3 2 1 4 5 6]);
% end


Y = Date(:,3);
M = Date(:,2);
D = Date(:,1);

[Lines, Rows] = size(Date);
switch Rows
 case 3
    F = zeros(Lines,1);
 case 4
    F = Date(:,4);
 case 6
    F = celestial.coo.convertdms(Date(:,4:6),'H','f');
 otherwise
    error('Illegal number of column in Date matrix');
end

B  = zeros(Lines,1);
Im3 = find(M<3);
Y(Im3) = Y(Im3) - 1;
M(Im3) = M(Im3) + 12;

Iy = find(Y>1582 | (Y==1582 & M>10) | (Y==1582 & M==10 & D>=15));
A = floor(Y(Iy).*0.01);
B(Iy) = 2 - A + floor(0.25.*A);
JD = floor(365.25.*(Y + 4716)) + floor(30.6001.*(M + 1)) + D + B - 1524.5 + F;

if (nargin>1)
   switch lower(Output)
    case 'jd'
       % do nothing
    case 'mjd'
       JD = JD - 2400000.5;
    otherwise
       error('Unknown Output option');
   end
end
