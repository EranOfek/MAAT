function JD=date2jd(Date,Output)
    % Convert Julian/Gregorian date to Julian Day
    % Package: @Time
    % Description: Convert Julian/Gregorian date to Julian Day.
    % Input  : - Gregorian of Julian date in one of the following
    %            formats:
    %            [Y, M, D, Frac]
    %            [Y, M, D]
    %            [Y] - first day of year
    %            [Y, M] - first day of month
    %            [Y, M, D, H M S]
    %            [Day, Month, Year, Day_Fraction]
    %            or [Day, Month, Year, Hour, Min, Sec]
    %            Alternatively this can be a string or a cell array of strings
    %            in which each string contains the date in the format:
    %            'yyyy-mm-ddTHH:MM:SS' (e.g., '2010-08-11T15:01:56')
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
    % Example: celestial.time.date2jd([1 1 2000 10 30 0]);
    %          celestial.time.date2jd([1 1 2000; 2 2 2000]);
    %          celestial.time.date2jd({'20101012T101010.111'});
    %          celetial.time.date2jd;  % JD of now (UTC)
    % Reliable: 1
    %--------------------------------------------------------------------------

    if nargin<2
        Output = 'jd';
        if nargin<1
            Date = [];
        end
    end

    if isempty(Date)
        Date = datevec(datetime('now', 'TimeZone', 'UTC'));
    end
    if ischar(Date) || iscellstr(Date)
        Date = celestial.time.str2date(Date);
    end

    if isnumeric(Date)
        [Nrow, Ncol] = size(Date);
        switch Ncol
            case 1
                % assume first day of year
                Date = [Date, ones(Nrow,2)];
                Frac = zeros(Nrow,1);
            case 2
                % assume first day of month
                Date = [Date, ones(Nrow,1)];
                Frac = zeros(Nrow,1);
            case 3
                % assume 0 hour
                % Date already has 3 columns
                Frac = zeros(Nrow,1);
            case 4
                % assume fraction of day in 4th column
                Frac = Date(:,4);
                Date = Date(:,1:3);

            case 6
                % assume Y M D H M S
                % convert HMS to Frac
                Frac = (Date(:,4).*3600 + Date(:,5).*60 + Date(:,6))./Time.SecInDay;
                Date = Date(:,1:3);
            otherwise
                error('Unknown Date option');
        end
    end
    Nrow = size(Date,1);

    Y = Date(:,1);
    M = Date(:,2);
    D = Date(:,3);

    B   = zeros(Nrow,1);
    Fm3 = M<3; 
    Y(Fm3) = Y(Fm3) - 1;
    M(Fm3) = M(Fm3) + 12;

    Fy = (Y>1582 | (Y==1582 & M>10) | (Y==1582 & M==10 & D>=15));
    A = floor(Y(Fy).*0.01);
    B(Fy) = 2 - A + floor(0.25.*A);
    JD = floor(365.25.*(Y + 4716)) + floor(30.6001.*(M + 1)) + D + B - 1524.5 + Frac;

    switch lower(Output)
        case 'jd'
            % do nothing
        case 'mjd'
            JD = JD - 2400000.5;
        otherwise
            error('Unknown Output option');
    end

end % jd function
