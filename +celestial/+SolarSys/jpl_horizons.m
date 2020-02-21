function [Cat]=jpl_horizons(varargin)
% Get JPL horizons ephemeris for a solar system body.
% Package: celestial.SolarSys
% Description: Get JPL horizons ephemeris for a solar system body.
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cat]=celestial.SolarSys.jpl_horizons;
%          [Cat]=celestial.SolarSys.jpl_horizons('ObjectInd','9804','StartJD',celestial.time.julday([14 6 2018]),'StopJD',  celestial.time.julday([20 6 2018]));
% Reliable: 2
%--------------------------------------------------------------------------

CatField      = AstCat.CatField;
ColCellField  = AstCat.ColCellField;
ColUnitsField = AstCat.ColUnitsField;


DefV.ObjectInd           = '499;';   % semicolumn telss horizon its a small body - e.g., '499;'
DefV.StartJD             = 2451545.5;
DefV.StopJD              = 2451555.5;
DefV.DateFormat          = 'JD';
DefV.StepSize            = 1;
DefV.StepSizeUnits       = 'd';
DefV.OutputColumns       = '1,9,10,13,19,20,23,24';  % https://ssd.jpl.nasa.gov/horizons.cgi?s_tset=1#top
DefV.OutCoo              = 'rad';
DefV.CENTER              = '500'; %'@sun'; %code for observer location. Earth -  '500', GAIA - '500@-139479', '675' - Palomar
DefV.WebOptions          = weboptions;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


BaseURL = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?';

I = 0;
I = I + 1;
Str(I).command = 'batch';
Str(I).value   = 1;

I = I + 1;
Str(I).command = 'COMMAND';
Str(I).value   = InPar.ObjectInd;

I = I + 1;
Str(I).command = 'CENTER';
Str(I).value   = InPar.CENTER;


I = I + 1;
Str(I).command = 'MAKE_EPHEM';
Str(I).value   = 'YES';

I = I + 1;
Str(I).command = 'TABLE_TYPE';
Str(I).value   = 'OBSERVER';

I = I + 1;
Str(I).command = 'START_TIME';
Str(I).value   = convert.time(InPar.StartJD,InPar.DateFormat,'StrDateO');
Str(I).value   = Str(I).value{1};

I = I + 1;
Str(I).command = 'STOP_TIME';
Str(I).value   = convert.time(InPar.StopJD,InPar.DateFormat,'StrDateO');
Str(I).value   = Str(I).value{1};

I = I + 1;
Str(I).command = 'STEP_SIZE';
Str(I).value   = sprintf('%d%%20%s',InPar.StepSize,InPar.StepSizeUnits);

I = I + 1;
Str(I).command = 'QUANTITIES';
Str(I).value   = InPar.OutputColumns;

I = I + 1;
Str(I).command = 'CSV_FORMAT';
Str(I).value   = 'YES';



Nstr       = numel(Str);
AllCommand = '';
for Istr=1:1:Nstr
    if (isnumeric(Str(Istr).value))
        Value = sprintf('%d',Str(Istr).value);
    else
        Value = sprintf('''%s''',Str(Istr).value);
    end
    
    if (Istr==1)
        AllCommand = sprintf('%s=%s',Str(Istr).command,Value);
    else
        AllCommand = sprintf('%s&%s=%s',AllCommand,Str(Istr).command,Value);
    end
end

UrlCommand = sprintf('%s%s',BaseURL,AllCommand);
Data = webread(UrlCommand,InPar.WebOptions);

% read header
Lines = regexp(Data,'\n','split');
Iline = find(~Util.cell.isempty_cell(strfind(Lines,'Date__')));
ColCell = regexp(Lines(Iline),',','split');
if isempty(ColCell)
    Cat = AstCat;
   return  
end
ColCell = Util.string.spacedel(ColCell{1});
Ncol    = numel(ColCell);

% read table
Tmp    = regexp(Data,'\$\$SOE(?<data>.+)\$\$EOE','names');
Table  = Tmp.data;
Format = Util.string.str_duplicate('%s ',Ncol);
Format = sprintf('%s\n',Format);
C = textscan(Table,Format,'Delimiter',',');

% remove empty columns
Flag = cellfun(@isempty,ColCell);
ColCell = ColCell(~Flag);
C       = C(~Flag);
Ncol    = numel(ColCell);

switch lower(InPar.OutCoo)
    case 'deg'
        OutCooUnits = 'd';
    case {'rad','radian'}
        OutCooUnits = 'r';
    otherwise
        error('Unkinown OutCoo option');
end

% reformat columns
for Icol=1:1:Ncol
    
    switch ColCell{Icol}
        case 'Date__(UT)__HR:MN'
            JD = datenum(C{Icol},'yyyy-mmm-dd HH:MM')-datenum('2000-01-01') + celestial.time.julday([1 1 2000]);
            
            C{Icol}        = JD;
            ColCell{Icol}  = 'JD';
            ColUnits{Icol} = 'JD';
            
        case {'R.A._(ICRF/J2000.0)','R.A._(ICRF)'}
            C{Icol} = celestial.coo.convertdms(C{Icol},'SHb',OutCooUnits);
            ColCell{Icol}  = 'RA';
            ColUnits{Icol} = 'OutCooUNits';
        
        case {'DEC_(ICRF/J2000.0)','DEC__(ICRF)'}
            C{Icol} = celestial.coo.convertdms(C{Icol},'SDb',OutCooUnits);
            ColCell{Icol}  = 'Dec';
            ColUnits{Icol} = 'OutCooUNits';
        case 'APmag'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol} = 'APmag';
            ColUnits{Icol} = 'mag';
        case 'S-brt'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'SurfMag';
            ColUnits{Icol} = 'mag/arcsec^2';
            
        case 'Illu%'
            %Fraction of target circular disk illuminated by Sun (phase), as seen by
            % observer.  Units: PERCENT
            C{Icol} = str2double(C{Icol})./100;
            ColCell{Icol}  = 'IllumFrac';
            ColUnits{Icol} = '';
            
        case 'Ang-diam'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'AngDiameter';
            ColUnits{Icol} = 'arcsec';
            
        case 'r'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'r';
            ColUnits{Icol} = 'au';
            
        case 'rdot'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'rdot';
            ColUnits{Icol} = 'km/s';
            
        case 'delta'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'Delta';
            ColUnits{Icol} = 'au';
        case 'deldot'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'Deltadot';
            ColUnits{Icol} = 'km/s';
        case 'S-O-T'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'SubObsTargetAng';
            ColUnits{Icol} = 'deg';
        case '/r'
            % 
            % -1 -- /T indicates target TRAILS Sun (evening sky; rises and sets AFTER Sun)
            % +1 -- /L indicates target LEADS Sun  (morning sky; rises and sets BEFORE Sun)

            FlagT = strcmp(C{Icol},'/T');
            Tmp   = zeros(numel(C{Icol}),1);
            Tmp(FlagT)  = -1;
            Tmp(~FlagT) = 1;
            
            C{Icol} = Tmp;
            ColCell{Icol}  = 'TrailLeadSun';
            ColUnits{Icol} = '';
        case 'S-T-O'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'SubTargetObsAng';
            ColUnits{Icol} = 'deg';
        case '1-way_LT'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'LightTime_1w';
            ColUnits{Icol} = 'minute';
        case 'ObsEcLon'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'ObsEcLon';
            ColUnits{Icol} = 'deg';        
        case 'ObsEcLat'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'ObsEcLat';
            ColUnits{Icol} = 'deg';      
        case 'T-mag'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'T_mag';
            ColUnits{Icol} = 'mag';      
	    case 'N-mag'
            C{Icol} = str2double(C{Icol});
            ColCell{Icol}  = 'N_mag';
            ColUnits{Icol} = 'mag';      
        otherwise
            ColCell{Icol}
            error('Unknwon column name option');
    end
end

Cat = AstCat;
Cat.(CatField) = [C{:}];
Cat.(ColCellField)  = ColCell;
Cat.(ColUnitsField) = ColUnits;
Cat = colcell2col(Cat);
