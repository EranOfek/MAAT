function [OutData,OutRad]=convertdms1(InData,InUnits,OutUnits)
%--------------------------------------------------------------------------
% convertdms1 function                                               ephem
% Description: 
% Input  : - 
% Output : - 
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------
import celestial.coo.*
import Util.cell.*

InvRAD = pi./180;
RAD    = 180./pi;

IsCell = false;
if (iscell(InData)),
    IsCell = true;
else
    if (ischar(InData)),
        InData = {InData};
    end
end
        
% treat gD and gH options
switch InUnits
    case 'gD'
        if (iscellstr(InData)),
            Tmp    = InData{1};
            IsNum = false;
        elseif (ischar(InData)),
            Tmp    = InData;
            IsNum = false;
        else
            Tmp    = InData(1,:); % InData is numeric
            IsNum = true;
        end
        
        if (IsNum),
            if (length(Tmp)==1),
                InUnits = 'r';
            elseif (length(Tmp)==2),
                InUnits = 'deg';
            elseif (length(Tmp)==3),
                InUnits = 'Dm';
            elseif (length(Tmp)==4),             
                InUnits = 'D';
            else
                error('Unknown InData format while gD option');
            end
            
        else
            if (length(strfind(Tmp,':'))==2 || ...
                length(strfind(Tmp,'s'))==1 || ...
                length(strfind(Tmp,' '))==2),
                InUnits = 'SD';
            elseif (length(strfind(Tmp,':'))==1 || ...
                    (length(strfind(Tmp,'m'))==1 && isempty(strfind(Tmp,'s')))|| ...
                    length(strfind(Tmp,' '))==1),
                % only minutes
                InUnits = 'SDm';
            else
                % assume no spaces
                InUnits = 'SDn';
            end
        end
    case 'gH'
        if (iscellstr(InData)),
            Tmp    = InData{1};
            IsNum = false;
        elseif (ischar(InData)),
            Tmp    = InData;
            IsNum = false;
        else
            Tmp    = InData(1,:); % InData is numeric
            IsNum = true;
        end
        
        if (IsNum),
            if (length(Tmp)==1),
                InUnits = 'r';
            elseif (length(Tmp)==2),
                InUnits = 'Hm';
            elseif (length(Tmp)==3),
                InUnits = 'H';
            else
                error('Unknown InData format while gD option');
            end
        
        else
            if (length(strfind(Tmp,':'))==2 || ...
                length(strfind(Tmp,'s'))==1 || ...
                length(strfind(Tmp,' '))==2),
                InUnits = 'SH';
            
            elseif (length(strfind(Tmp,':'))==1 || ...
                    (length(strfind(Tmp,'m'))==1 && isempty(strfind(Tmp,'s')))|| ...
                    length(strfind(Tmp,' '))==1),
                % only minutes
                InUnits = 'SHm';
            else
                % assume no spaces
                InUnits = 'SHn';
            end
        end   
          
    otherwise
                % do nothing
end
            
          
     InUnits       
% conver InData to radians            
switch InUnits
    case {'r','R'}
        % already in radians - do nothing
        OutRad = InData;
    case 'f'
        % convert fraction of day to radians in the 0..2pi range
        OutRad = InData.*2.*pi;
    case {'d','deg'}
        % convert degrees to radians
        OutRad = InData.*InvRAD;
    case 'D'
        % convert [Sign D M S] to radians
        OutRad = InData(:,1).*(InData(:,2) + InData(:,3)./60 + InData(:,4)./3600).*InvRAD;
    case {'DM','Dm','Dmin'}
        % convert [Sign D M] to radians
        OutRad = InData(:,1).*(InData(:,2) + InData(:,3)./60).*InvRAD;
    case 'h'
        % convert hours to radians
        OutRad = InData.*15.*InvRAD;
    case 'H'
        % convert [H M S] to radians
        OutRad = (InData(:,1) + InData(:,2)./60 + InData(:,3)./3600).*15.*InvRAD;
    case {'HM','Hmin'}
        % convert [H M] to radians
        OutRad = (InData(:,1) + InData(:,2)./60).*15.*InvRAD;
    case {'SD','SDh','SDb','SH','SHh','SHb'}
        % convert sexagesimal string to radians
        Nin = numel(InData);
        for Iin=1:1:Nin,
            DMS(Iin) = regexp(InData{Iin},'(?<Sign>[\s\+-]?)(?<D>\d\d)[:hd\s](?<M>\d\d)[:m\s](?<S>\d\d\.*\d*)','names');
        end
        %Fp = strcmp({DMS.Sign},'+') | strcmp({DMS.Sign},'') | strcmp({DMS.Sign},' ');
        Fm = strcmp({DMS.Sign},'-');
        D  = str2double({DMS.D});
        M  = str2double({DMS.M});
        S  = str2double({DMS.S});
        Sign = ones(size(D));
        Sign(Fm) = -1;
        OutRad = [Sign.*(D + M./60 + S./3600).*InvRAD].';
 
        if (strcmp(InUnits(2),'H')),
            % Hours
            OutRad = OutRad.*15;
        end
        
            
    case {'SDn','SHn'}
        % convert sexagesimal string without spaces to radians
        Nin = numel(InData);
        for Iin=1:1:Nin,
            DMS(Iin) = regexp(InData{Iin},'(?<Sign>[\s\+-]?)(?<D>\d\d)(?<M>\d\d)(?<S>\d\d\.*\d*)','names');
        end
        %Fp = strcmp({DMS.Sign},'+') | strcmp({DMS.Sign},'') | strcmp({DMS.Sign},' ');
        Fm = strcmp({DMS.Sign},'-');
        D  = str2double({DMS.D});
        M  = str2double({DMS.M});
        S  = str2double({DMS.S});
        Sign = ones(size(D));
        Sign(Fm) = -1;
        OutRad = [Sign.*(D + M./60 + S./3600).*InvRAD].';
        
        if (strcmp(InUnits(2),'H')),
            % Hours
            OutRad = OutRad.*15;
        end
    case {'SDm','SHm'}
        % convert sexagesimal strings without seconds to radians
        Nin = numel(InData);
        for Iin=1:1:Nin,
            DMS(Iin) = regexp(InData{Iin},'(?<Sign>[\s\+-]?)(?<D>\d\d)[:hd\s](?<M>\d\d\.*\d*)','names');
        end
        %Fp = strcmp({DMS.Sign},'+') | strcmp({DMS.Sign},'') | strcmp({DMS.Sign},' ');
        Fm = strcmp({DMS.Sign},'-');
        D  = str2double({DMS.D});
        M  = str2double({DMS.M});
        Sign = ones(size(D));
        Sign(Fm) = -1;
        OutRad = [Sign.*(D + M./60).*InvRAD].';
        
        if (strcmp(InUnits(2),'H')),
            % Hours
            OutRad = OutRad.*15;
        end
    otherwise
        error('unknown InUnits option');
end




switch OutUnits
    case {'r','R'}
        % convert to radians
        OutData = OutRad;
    case 'f'
        % convert to fraction of day
        OutData = OutRad./(2.*pi);
    case {'d','deg'}
        % convert to deg
        OutData = OutRad.*RAD;
    case {'D','SD','SDh','SDb','SDn'}
        % convert to [Sign D M S]
        OutData = zeros(length(OutRad),4);
        Tmp = abs(OutRad.*RAD);
        OutData(:,1:2) = [sign(OutRad), floor(Tmp)];
        OutData(:,3) = (Tmp-OutData(:,2)).*60;
        OutData(:,4) = (OutData(:,3) - floor(OutData(:,3))).*60;
        OutData(:,3) = floor(OutData(:,3));
        
        switch OutUnits
            case 'SD'
                ...OutData = sprintf2cell('%02d:%02d:%06.3f',OutData);
            case 'SDh'
                ...OutData = sprintf2cell('%02dh%02dm%06.3fs',OutData);
            case 'SDb'
                ...OutData = sprintf2cell('%02d %02d %06.3f',OutData);
            case 'SDn'
                ...OutData = sprintf2cell('%02d%02d%06.3f',OutData);
            otherwise
                % do nothing
        end
    case {'DM','Dm','Dmin'}
        % convert to [Sign D M]
        OutData = zeros(length(OutRad),4);
        Tmp = abs(OutRad.*RAD);
        OutData(:,1:2) = [sign(OutRad), floor(Tmp)];
        OutData(:,3) = floor((Tmp-OutData(:,2)).*60);
        
    case 'h'
        % convert to hours
        OutData = OutRad.*RAD./15;
        
    case {'H','SH','SHh','SHb','SHn'}
        % convert to [H M S]
        OutData = zeros(length(OutRad),3);
        Tmp = abs(OutRad.*RAD./15);
        OutData(:,1) = floor(Tmp);
        OutData(:,2) = (Tmp-OutData(:,1)).*60;
        OutData(:,3) = (OutData(:,2) - floor(OutData(:,2))).*60;
        OutData(:,2) = floor(OutData(:,2));
        
        switch OutUnits
            case 'SH'
                OutData = sprintf2cell('%02d:%02d:%06.3f',OutData);
            case 'SHh'
                OutData = sprintf2cell('%02dh%02dm%06.3fs',OutData);
            case 'SHb'
                OutData = sprintf2cell('%02d %02d %06.3f',OutData);
            case 'SHn'
                OutData = sprintf2cell('%02d%02d%06.3f',OutData);
            otherwise
                % do nothing
        end
    case {'HM','Hmin','SHm'}
        % convert to [H M]
        OutData = zeros(length(OutRad),2);
        Tmp = abs(OutRad.*RAD./15);
        OutData(:,1) = floor(Tmp);
        OutData(:,2) = (Tmp-OutData(:,1)).*60;
        
        switch OutUnits
            case 'SHm'
                OutData = sprintf2cell('%02d:%06.3f',OutData);
            otherwise
                % do nothing
        end
        
    case {'SD','SDh','SDb','SH','SHh','SHb'}
        
    case {'SDn','SHn'}
       
    case {'SDm','SHm'}
      
    otherwise
        error('unknown InUnits option');
end
