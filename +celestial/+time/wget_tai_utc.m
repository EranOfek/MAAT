function TAI_UTC=wget_tai_utc(Read)
% Get TAI-UTC from file or IERS website
% Package: celestial.time
% Description: Get the table of historical TAI-UTC time differences 
%              (leap second) from the IERS web site.
% Input  : - 'get' - get the latest TAI-UTC data from the IERS website and
%                    update the local version.
%            'use' - use local version of the TAI-UTC data (default).
% Output : - A structure containing the TAI-UTC data. The structure
%            contains the following fields:
%            .Cat  - The catalog.
%            .Col  - A structure describing the catalog columns.
%            .UnitsCell - A cell array of the units of each column.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: TAI_UTC=celestial.time.wget_tai_utc;
% Reliable: 2
%--------------------------------------------------------------------------


Def.Read = 'use';
if (nargin==0)
   Read = Def.Read;
end

Dir     = Util.files.which_dir('celestial.time.wget_tai_utc');
Dir     = sprintf('%s%s..%s..%s..%sdata%sSolarSystem%s',Dir,filesep,filesep,filesep,filesep,filesep,filesep);

switch lower(Read)
    case 'get'
        % source
        % http://maia.usno.navy.mil/ser7/tai-utc.dat


        URL = 'http://maia.usno.navy.mil/ser7/tai-utc.dat';
        Str = urlread(URL);

        %Res=regexp(Str,'=JD (?<JD>\d+\.\d)','names')
        %Res=regexp(Str,'TAI-UTC=\s+(?<TU>\d+\.\d+)','names')
        %Res=regexp(Str,'MJD - (?<T0>\d+\.\d*)','names')
        %Res=regexp(Str,') X (?<S>\d+\.\d*)','names')

        Res = regexp(Str,'=JD (?<JD>\d+\.\d)\s+TAI-UTC=\s+(?<TU>\d+\.\d+)\s+S \+ \(MJD - (?<T0>\d+\.\d*)\) X (?<S>\d+\.\d*)','names');
        JD         = str2double({Res.JD});
        DelTAU_UTC = str2double({Res.TU});
        T0         = str2double({Res.T0});
        S          = str2double({Res.S});

        TAI_UTC.Cat = [JD.', DelTAU_UTC.', T0.', S.'];
        TAI_UTC.Col.JD         = 1;
        TAI_UTC.Col.DelTAU_UTC = 2;
        TAI_UTC.Col.T0         = 3;
        TAI_UTC.Col.S          = 4;
        TAI_UTC.UnitsCell      = {'day','s','day','s/day'};
        TAI_UTC.Fun            = @tai_utc;
        
        % save data
        save(sprintf('%sTAI_UTC.mat',Dir),'TAI_UTC');
    case 'use'
        TAI_UTC = Util.IO.load2(sprintf('%sTAI_UTC.mat',Dir));
    otherwise
        error('Unknown Read option');
end

       
