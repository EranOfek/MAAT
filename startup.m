%function startup(AddRem,ComputerPath)
% startup script
% Package: null
% Description: matlab startup function and set the matlab default behaviour.
%              This is including randomizing the seed of the matlab random
%              number generator and setting the figures font defaults.
% Input  : - Add or remove the astronomy & astrophysics matlab packages to
%            the matlab path {1 - for add|0 - for remove}. Default is 1 (add).
%          - Root Directory (check code for default).
% Output : null
% Tested : Matlab R2011a
%     By : Eran O. Ofek                    Sep 1996
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: startup('n'); % remove from the matlab path variable.
%--------------------------------------------------------------------------


% Assign variables to the matlab workspace
assignin('base','RAD',180./pi);    % Radian

% set the plot AxesFontSize and AxesFontName default
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultAxesFontName','times');

format short g

% randomizing the seed of the matlab random number generator
%rand('state',sum(100*clock));
rng('shuffle');



% Set the location of the functions and catalogs directories
%--- EDIT these lines ---

if (ismac || isunix)
    % mac / unix
    HomeDir = getenv('HOME');
    
else
    % windows
    HomeDir = getenv('HOMEPATH');
end

if (isempty(HomeDir))
    error('Can not find home directory environment variable - edit the startup.m file accordingly');
end

BaseDir = 'matlab';

ComputerPath = sprintf('%s%s%s%s',HomeDir,filesep,BaseDir,filesep);  % e.g., '/home/eran/matlab'
%---

%RunStartup = true;  % set to false if you want to avoid running this file


%
ListDir = {'MAAT',...
           'MAAT/FRT',...
           'MAAT/Catalogue',...
           'data',...
           'data/BATSE_LC',...
           'data/ELP',...
           'data/SolarSystem',...
           'data/VSOPE87'};


% ListDirCats= {'2MASS',...
%               '2MASSxsc',...
% 	      'AAVSO_VSX',...
%               'AKARI',...
%               'APASS',...
% 	      'Cosmos',...
% 	      'CRTS_per_var',...
%               'DECaLS/DR5',...
%               'FIRST',...
%               'GAIA/DR1',...
%               'GAIA/DR2',...
%               'GALEX/DR6Plus7',...
% 	      'GLIMPSE',...
%               'HST/HSCv2',...
%               'IPHAS/DR2',...
% 	      'LAMOST/DR4',...
%               'NED/20180502',...
%               'NVSS',...
% 	      'PGC',...
%               'PS1',...
%               'PTFpc',...
%               'ROSATfsc',...
% 	      'SDSS/DR10',...
%               'SDSS/DR14offset',...
% 	      'Simbad_PM200',...
% 	      'SkyMapper',...
% 	      'SpecSDSS/DR14',...
% 	      'Spitzer/SAGE',...
% 	      'Spitzer/IRACgc',...
%               'SWIREz',...
%               'UCAC4',...
%               'UKIDSS/DR10',...
% 	      'URAT1',...
%               'USNOB1',...
%               'VISTA/Viking/DR2',...
%               'VST/ATLAS/DR3',...
%               'VST/KiDS/DR3',...
%               'WISE',...
%               'XMM'};

ListDirCats= {'2MASS',...
              '2MASSxsc',...
	      'AAVSO_VSX',...
              'AKARI',...
              'APASS',...
	      'Cosmos',...
	      'CRTS_per_var',...
              'DECaLS/DR5',...
              'FIRST',...
              'GAIA/DR1',...
              'GAIA/DR2',...
              'GALEX/DR6Plus7',...
	      'GLIMPSE',...
              'HST/HSCv2',...
              'IPHAS/DR2',...
	      'LAMOST/DR4',...
              'NED/20180502',...
              'NVSS',...
	      'PGC',...
              'PS1',...
              'PTFpc',...
              'ROSATfsc',...
	      'SDSS/DR10',...
              'SDSS/DR14offset',...
	      'Simbad_PM200',...
	      'SkyMapper',...
	      'SpecSDSS/DR14',...
	      'Spitzer/SAGE',...
	      'Spitzer/IRACgc',...
              'SWIREz',...
              'UCAC4',...
              'UKIDSS/DR10',...
	      'URAT1',...
          'unWISE',...
              'USNOB1',...
              'VISTA/Viking/DR2',...
              'VST/ATLAS/DR3',...
              'VST/KiDS/DR3',...
              'WISE',...
              'XMM',...
              'ZTF/LCDR1',...
              'ZTF/SrcLCDR1'};

ListDir = strrep(ListDir,'/',filesep);
ListDirCats = strrep(ListDirCats,'/',filesep);

Add2path = ListDir;
for I=1:1:numel(ListDir)
    Add2path{I} = sprintf('%s%s',ComputerPath,Add2path{I});
end

CatsPath     = sprintf('/euler/catsHTM/');

Add2pathCats = ListDirCats;
for I=1:1:numel(ListDirCats)
    Add2pathCats{I} = sprintf('%s%s',CatsPath,Add2pathCats{I});
end



addpath(Add2path{:});
addpath(Add2pathCats{:});


clear HomeDir      
clear BaseDir
clear ComputerPath
clear CatsPath
clear ListDir
clear ListDirCats
clear Add2path
clear Add2pathCats
