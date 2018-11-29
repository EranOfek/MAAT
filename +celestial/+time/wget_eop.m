function EOP=wget_eop(Read)
% Read or get the Earth orientation parameters file from IERS
% Package: celestial.time
% Description: Get the table of historical and predicted Earth orientation
%              parameters (EOP) from the IERS web site.
% Input  : - 'get' - get the latest EOP data from the IERS website and
%                    update the local version.
%            'use' - use local version of the EOP data (default).
% Output : - A structure containing the EOP data. The structure contains
%            the following fields:
%            .Cat  - The catalog.
%            .Col  - A structure describing the catalog columns.
%            .UnitsCell - A cell array of the units of each column.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: EOP=celestial.time.wget_eop;
% Reliable: 2
%--------------------------------------------------------------------------

Def.Read = 'use';
if (nargin==0)
   Read = Def.Read;
end

Dir     = Util.files.which_dir('celestial.time.wget_eop');
Dir     = sprintf('%s%s..%s..%s..%sdata%sSolarSystem%s',Dir,filesep,filesep,filesep,filesep,filesep,filesep);

switch lower(Read)
    case 'get'
      %  URL = 'https://datacenter.iers.org/eop/-/somos/5Rgv/latestXL/9/finals2000A.all/csv';
        URL = 'https://datacenter.iers.org/data/csv/finals2000A.data.csv';
        Str=urlread(URL);
        Format = '%f %f %f %f %s %f %f %f %f %s %f %f %f %f %s %f %f %f %f %f %f %f %f %s %f %f %s %f %s %f %f %f %f\n';
        Cell = textscan(Str,Format,'Delimiter',';','HeaderLines',1);

        % documentation: http://hpiers.obspm.fr/eoppc/bul/bulb/explanatory.html
        % http://maia.usno.navy.mil
        % http://www.iers.org/nn_10968/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true

        % columns:
        %MJD;Year;Month;Day;Type;X;sigma_X;Y;sigma_Y;Type;UT1-UTC;sigma_UT1-UTC;LOD;sigma_LOD;Type;dPsi;sigma_dPsi;dEpsilon;sigma_dEpsilon;dX;sigma_dX;dY;sigma_dY;Type;bulB/X;bulB/Y;Type;bulB/UT-UTC;Type;bulB/dPsi;bulB/dEpsilon;bulB/dX;bulB/dY

        Cell = Cell([1:4,6:9,11:14,16:23,25:26,28,30:33]);
        EOP.Cat = cell2mat(Cell);
        Col.MJD     = 1;
        Col.Year    = 2;
        Col.Month   = 3;
        Col.Day     = 4;
        %Col.Type    = 5;
        Col.X       = 5;  % Celestial Ephemeris Pole (CEP) relative to the International Reference Pole (IRP) are defined as x and y
        Col.sigmaX  = 6;
        Col.Y       = 7;
        Col.sigmaY  = 8;
        %Col.Type    = 10;
        Col.UT1_UTC = 9;
        Col.sigmaUT1_UTC = 10;   % UT1-UTC
        Col.LOD     = 11;
        Col.sigmaLOD= 12;
        %Col.Type    = 15;
        Col.dPsi    = 13;    % offset relative to IAU 1980 Theory of Nutation
        Col.sigmadPsi = 14;
        Col.dEps    = 15;    % % offset relative to IAU 1980 Theory of Nutation
        Col.sigmadEps = 16;
        Col.dX      = 17;
        Col.sigmadX = 18;
        Col.dY      = 19;
        Col.sigmadY = 20;
        %Col.Type    = 24;
        Col.bulBx   = 21;
        Col.bulBy   = 22;
        %Col.Type    = 27;
        Col.bulB_UT_UTC = 23;
        %Col.Type    = 29;
        Col.bulBdPsi = 24;
        Col.bulBdEps = 25;
        Col.bulBdX   = 26;
        Col.bulBdY   = 27;

        EOP.Col = Col;
        EOP.UnitsCell = {'day','','','','mas','mas','mas','mas','s','s','ms','ms','mas','mas','mas','mas','mas','mas','mas','mas','mas','mas','mas','mas','mas','mas','mas'};

        % save data
        save(sprintf('%sEOP.mat',Dir),'EOP');
    case 'use'
        EOP = Util.IO.load2(sprintf('%sEOP.mat',Dir));
    otherwise
        error('Unknown Read option');
end

       
