function Data=get_orbit_files(Read)
%------------------------------------------------------------------------------
% get_orbit_files function                                           Catalogue
% Description: Get asteroids and comets orbital elements from JPL
%              and read into a matlab structure.
% Input  : - 'get' - get the latest orbital elements file from JPL
%            'use' - use local orbital elements file (default).
% Output : - Data structure containing the following fields:
%            .Cat    : 1 - num. asteroid; 2 - unnum asteroid; 3 - comet.
%            .Number : asteroid number
%            .Name   : Designation
%            .Epoch  : Epoch [JD]
%            .a      : semi major axis [AU]
%            .q      : perihelion distance [AU]
%            .e      : eccentricity
%            .i      : J2000 Inclination [deg]
%            .w      : J2000 Argument of perihelion [deg] 
%            .Node   : J2000 Longitude of ascending node [deg]
%            .MeanAnom: Mean anomaly at Epoch [deg]
%            .P      : Period [s]
%            .n      : Mean motion [deg/day]
%            .Tp     : Time of perihelion [JD]
%            .H      : Asteroid abs. mag. [mag]
%            .G      : Magnitude slope parameter
%            .Ref    : Reference
% Tested : Matlab 7,8
%     By : Eran O. Ofek                  November 2009
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Files source:	http://ssd.jpl.nasa.gov/?sb_elem
% Example: Data=get_orbit_files('get');
% Reliable: 1
%------------------------------------------------------------------------------
Files{1}='ELEMENTS.NUMBR';
Files{2}='ELEMENTS.UNNUM';
Files{3}='ELEMENTS.COMET';


Def.Read = 'use';
if (nargin==0),
   Read = Def.Read;
end

G = get_constant('G');
AU = get_constant('au');
SolarMass = get_constant('SolM');
PWD = pwd;
fprintf('cd to ../matlab/data/SolarSystem/ directory\n');

Dir     = which_dir('get_orbit_files.m');
Dir     = sprintf('%s/../../data/SolarSystem/',Dir)

cd(Dir);

switch lower(Read)
 case 'get'
    % get orbital elements from JPL:
    delete(Files{1});
    delete(Files{2});
    delete(Files{3});
    delete(sprintf('%s_1',Files{1}));
    delete(sprintf('%s_1',Files{2}));
    delete(sprintf('%s_1',Files{3}));
    system('wget http://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR.gz');
    system('wget http://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM.gz');
    system('wget wget http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET');
    system('gzip -d ELEMENTS*.gz');
    system('cut -c 1-105 ELEMENTS.NUMBR > ELEMENTS.NUMBR_1');
    system('cut -c 1-93 ELEMENTS.UNNUM > ELEMENTS.UNNUM_1');
    system('cut -c 1-114 ELEMENTS.COMET > ELEMENTS.COMET_1');
 otherwise
    % do nothing
end

%--- read orbital elements of numbered asteroids ---
FID=fopen('ELEMENTS.NUMBR_1','r');
%C=textscan(FID,'%f %17c %f %f %f %f %f %f %f %f %f %10c','Headerlines',2);
C=textscan(FID,'%f %17c %f %f %f %f %f %f %f %f %f','Headerlines',2);
fclose(FID);
C{2}=cellstr(C{2});                                                       
%C{12}=cellstr(C{12});
N = length(C{1});
Data1.Cat    = ones(N,1).*1;
Data1.Number = C{1};
Data1.Name   = C{2};
Data1.Epoch  = C{3} + 2400000.5;
Data1.a      = C{4};
Data1.e      = C{5};
Data1.q      = Data1.a.*(1-Data1.e);
Data1.i      = C{6};
Data1.w      = C{7};
Data1.Node   = C{8};
Data1.MeanAnom = C{9};
Data1.H        = C{10};
Data1.G        = C{11};
Data1.P        = sqrt(4.*pi.^2.*(AU.*Data1.a).^3./(G.*SolarMass));   % [s]
Data1.n        = 360./(Data1.P./86400);  % [deg/day]
Data1.Tp     = Data1.Epoch - Data1.MeanAnom./Data1.n;
%Data1.Ref = C{12};

%--- read orbital elements of unnumbered asteroids ---
FID=fopen('ELEMENTS.UNNUM_1','r');
%C=textscan(FID,'%11c %f %f %f %f %f %f %f %f %f %10c','Headerlines',2);
C=textscan(FID,'%11c %f %f %f %f %f %f %f %f %f','Headerlines',2);
fclose(FID);
C{1}=cellstr(C{1});                                                       
%C{11}=cellstr(C{11});
N = length(C{1});
Data2.Cat    = ones(N,1).*2;
Data2.Number = zeros(N,1).*NaN;
Data2.Name   = C{1};
Data2.Epoch  = C{2} + 2400000.5;
Data2.a      = C{3};
Data2.e      = C{4};
Data2.q      = Data2.a.*(1-Data2.e);
Data2.i      = C{5};
Data2.w      = C{6};
Data2.Node   = C{7};
Data2.MeanAnom = C{8};
Data2.H        = C{9};
Data2.G        = C{10};
Data2.P        = sqrt(4.*pi.^2.*(AU.*Data2.a).^3./(G.*SolarMass));   % [s]
Data2.n        = 360./(Data2.P./86400);  % [deg/day]
Data2.Tp     = Data2.Epoch - Data2.MeanAnom./Data2.n;
%Data2.Ref = C{11};


%--- read orbital elements of comets ---
FID=fopen('ELEMENTS.COMET_1','r');
%C=textscan(FID,'%38c %f %f %f %f %f %f %f %12c','Headerlines',2);
C=textscan(FID,'%45c %f %f %f %f %f %f %f','Headerlines',2);
fclose(FID);
C{1}=cellstr(C{1});                                                       
%C{9}=cellstr(C{9});
N = length(C{1});
Data3.Cat    = ones(N,1).*3;
Data3.Number = zeros(N,1).*NaN;
Data3.MeanAnom = zeros(N,1).*NaN;
Data3.H        = zeros(N,1).*NaN;
Data3.G        = zeros(N,1).*NaN;
Data3.Name   = C{1};
Data3.Epoch  = C{2} + 2400000.5;
Data3.q      = C{3};
Data3.e      = C{4};
Data3.a      = Data3.q./(1-Data3.e);
Data3.i      = C{5};
Data3.w      = C{6};
Data3.Node   = C{7};
Frac  = abs(C{8})-floor(abs(C{8}));
Day   = mod(abs(fix(C{8})),100);
Month = mod((abs(fix(C{8}))-Day)./100,100);
Year  = ((abs(fix(C{8}))-Day)./100-Month)./100;
Data3.Tp     = julday([Day, Month, Year, Frac]);
%Data3.Ref = C{9};
Data3.P        = sqrt(4.*pi.^2.*(AU.*Data3.a).^3./(G.*SolarMass));
Data3.n        = 360./(Data3.P./86400);  % [deg/day]

% go back to original directory
cd(PWD);


Data = structcon(Data1,Data2,1);
Data = structcon(Data,Data3,1);
