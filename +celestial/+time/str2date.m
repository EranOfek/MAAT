function Date=str2date(String)
    % Convert a date string (usinf datevec) to date vector
    % Input : - Date string (or cell of strings) that may have one of the following formats:
    %           yyyymmddTHHMMSS.FFF
    %           yyyy:mm:dd HH:MM:SS.FFF
    %           yyyy-mm-dd HH:MM:SS.FFF
    % Output - A date vector [y m d, H M S]
    % Example: Date=celestial.time.str2date({'2020:02:01 10:10:23.1123'})

    Failed = true;
    if Failed
        try
            Date = datevec(String,'yyyymmddTHHMMSS.FFF');
            Failed = false;
        end
    end    
    if Failed
        try
            Date = datevec(String,'yyyy:mm:dd HH:MM:SS.FFF');
            Failed = false;
        end
    end   
    if Failed
        try
            Date = datevec(String,'yyyy-mm-dd HH:MM:SS.FFF');
            Failed = false;
        end
    end   
    if Failed
        Date = [];
    end

end