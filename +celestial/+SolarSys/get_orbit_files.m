function Data=get_orbit_files(Read,OutType)
% Get asteroids/comets orbital elements from JPL, save locally and read.
% Package: celestial.SolarSys
% Description: Get asteroids and comets orbital elements from JPL
%              and read into a matlab structure.
% Input  : - 'wget','get' - get the latest orbital elements file from JPL
%            'load','use' - use local orbital elements file (default).
%          - Output type: 'struct'|'orbitalel'. Default is 'orbitalel'
% Output : - Data structure containing the following fields:
%            .Cat    : 1 - num. asteroid; 2 - unnum asteroid; 3 - comet.
%            .Number : asteroid number
%            .Name   : Designation
%            .Epoch  : Epoch [JD]
%            .a      : semi major axis [AU]
%            .q      : perihelion distance [AU]
%            .e      : eccentricity
%            .i      : J2000 Inclination [rad]
%            .w      : J2000 Argument of perihelion [rad] 
%            .Om     : J2000 Longitude of ascending node [rad]
%            .M      : Mean anomaly at Epoch [rad]
%            .P      : Period [s]
%            .n      : Mean motion [deg/day]
%            .Tp     : Time of perihelion [JD]
%            .H      : Asteroid abs. mag. [mag]
%            .G      : Magnitude slope parameter
%            .Ref    : Reference
% Tested : Matlab 7,8
%     By : Eran O. Ofek                    Nov 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Files source:	http://ssd.jpl.nasa.gov/?sb_elem
% Example: Data=celestial.SolarSys.get_orbit_files('get');
% Reliable: 1
%--------------------------------------------------------------------------
Files{1}='ELEMENTS.NUMBR';
Files{2}='ELEMENTS.UNNUM';
Files{3}='ELEMENTS.COMET';

if (nargin<2)
    OutType = 'orbitalel';
end

Def.Read = 'use';
if (nargin==0)
   Read = Def.Read;
end

RAD = 180./pi;
G   = constant.G; %get_constant('G');
AU  = constant.au; %get_constant('au');
SolarMass = constant.SunM; %get_constant('SolM');
PWD = pwd;
fprintf('cd to ../matlab/data/SolarSystem/ directory\n');

%Dir     = which_dir('get_orbit_files.m');
%Dir     = sprintf('%s/../../data/SolarSystem/',Dir);
Dir     = sprintf('~/matlab/data/SolarSystem/');

cd(Dir);

switch lower(Read)
 case {'wget','get'}
    % get orbital elements from JPL:
    delete(Files{1});
    delete(Files{2});
    delete(Files{3});
    delete(sprintf('%s_1',Files{1}));
    delete(sprintf('%s_1',Files{2}));
    delete(sprintf('%s_1',Files{3}));
    www.pwget('http://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR.gz','-nc');
    www.pwget('http://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM.gz','-nc');
    www.pwget('http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET','-nc');
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

Data1.ObjType= ones(N,1).*1;
Data1.Number = C{1};
Data1.Name   = C{2};
Data1.Epoch  = C{3} + 2400000.5;
a      = C{4};
Data1.e      = C{5};
%Data1.q      = Data1.a.*(1-Data1.e);
Data1.q      = C{4}.*(1-Data1.e);
Data1.i      = C{6}./RAD;
Data1.w      = C{7}./RAD;
Data1.Om     = C{8}./RAD;
Data1.M      = C{9}./RAD;
Data1.H        = C{10};
Data1.G        = C{11};
P        = sqrt(4.*pi.^2.*(AU.*a).^3./(G.*SolarMass));   % [s]
n        = 2.*pi./(P./86400);  % [rad/day]
Data1.T      = Data1.Epoch - Data1.M./n;
%Data1.Ref = C{12};

%--- read orbital elements of unnumbered asteroids ---
FID=fopen('ELEMENTS.UNNUM_1','r');
%C=textscan(FID,'%11c %f %f %f %f %f %f %f %f %f %10c','Headerlines',2);
C=textscan(FID,'%11c %f %f %f %f %f %f %f %f %f','Headerlines',2);
fclose(FID);
C{1}=cellstr(C{1});                                                       
%C{11}=cellstr(C{11});
N = length(C{1});

Data2.ObjType= ones(N,1).*2;
Data2.Number = zeros(N,1).*NaN;
Data2.Name   = C{1};
Data2.Epoch  = C{2} + 2400000.5;
a      = C{3};
Data2.e      = C{4};
%Data2.q      = Data2.a.*(1-Data2.e);
Data2.q      = C{3}.*(1-Data2.e);
Data2.i      = C{5}./RAD;
Data2.w      = C{6}./RAD;
Data2.Om     = C{7}./RAD;
Data2.M      = C{8}./RAD;
Data2.H        = C{9};
Data2.G        = C{10};
P        = sqrt(4.*pi.^2.*(AU.*a).^3./(G.*SolarMass));   % [s]Data2.n        = 360./(Data2.P./86400);  % [deg/day]
n        = 2.*pi./(P./86400);  % [rad/day]
Data2.T      = Data2.Epoch - Data2.M./n;
%Data2.Ref = C{11};


%--- read orbital elements of comets ---
FID=fopen('ELEMENTS.COMET_1','r');
%C=textscan(FID,'%38c %f %f %f %f %f %f %f %12c','Headerlines',2);
C=textscan(FID,'%40c %f %f %f %f %f %f %f','Headerlines',2);
fclose(FID);
C{1}=cellstr(C{1});                                                       
%C{9}=cellstr(C{9});
N = length(C{1});

Data3.ObjType= ones(N,1).*3;
Data3.Number = zeros(N,1).*NaN;
Data3.M      = zeros(N,1).*NaN;
Data3.H        = zeros(N,1).*NaN;
Data3.G        = zeros(N,1).*NaN;
Data3.Name   = C{1};
Data3.Epoch  = C{2} + 2400000.5;
Data3.q      = C{3};
Data3.e      = C{4};
a      = Data3.q./(1-Data3.e);
Data3.i      = C{5}./RAD;
Data3.w      = C{6}./RAD;
Data3.Om     = C{7}./RAD;
Frac  = abs(C{8})-floor(abs(C{8}));
Day   = mod(abs(fix(C{8})),100);
Month = mod((abs(fix(C{8}))-Day)./100,100);
Year  = ((abs(fix(C{8}))-Day)./100-Month)./100;
Data3.T      = celestial.time.julday([Day, Month, Year, Frac]);
%Data3.Ref = C{9};
P        = sqrt(4.*pi.^2.*(AU.*a).^3./(G.*SolarMass));
n        = 2.*pi./(P./86400);  % [rad/day]

% go back to original directory
cd(PWD);

Data = Util.struct.structcon(Data1,Data2,1);
Data = Util.struct.structcon(Data,Data3,1);

switch lower(OutType)
    case 'orbitalel'
        % convert struct to OrbitalEl
        Data = OrbitalEl.struct2orbitalel(Data);
end