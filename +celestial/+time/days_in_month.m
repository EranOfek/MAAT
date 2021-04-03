
function Nd=days_in_month(Year,Month)
    % Return the number of days in month
    % Package: @Time, adapted from +celestial.time
    % Description: Calculate the number of days in a given Gregorian or Julian
    %              month.
    % Input  : - Year (array)
    %          - Month (array, the same size as Year).
    % Output : - Number of days in month. 
    % Tested : Matlab 5.3
    %     By : Eran O. Ofek                    Jan 2003
    %    URL : http://weizmann.ac.il/home/eofek/matlab/
    % Example: Nd=celestial.time.days_in_month(2000,2);
    % Reliable: 2

    if numel(Year)~=numel(Month)
        error('Month and Year must have the same size');
    end
    MonthE = Month + 1;
    YearE  = Year;
    Flag = MonthE==13;
    MonthE(Flag)  = 1;
    YearE(Flag)    = Year(Flag)+1;

    Size = size(Month);
    N = numel(Year);
    JD1 = convert.date2jd([ones(N,1) Month(:)  Year(:)  zeros(N,1)]);
    JD2 = convert.date2jd([ones(N,1) MonthE(:) YearE(:) zeros(N,1)]);

    Nd = JD2 - JD1;

    Nd = reshape(Nd,Size);
end
