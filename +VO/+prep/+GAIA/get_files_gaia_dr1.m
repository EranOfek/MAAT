function get_files_gaia_dr1
% Get GAIA DR1 files from GAIA archive
% Package: VO.prep.GAIA
% Description: Get GAIA DR1 files from GAIA archive
% Input  : null (see internal parameters)
% Output : null
% Example: VO.prep.GAIA.get_files_gaia_dr1
% Reliable: 2

Type         = 'csv';   % 'csv'|'fits'|'votable'
BaseFileName = 'GaiaSource';
BaseURL      = 'http://cdn.gea.esac.esa.int/Gaia/gaia_source/';
CopyTo       = '/raid/eran/Catalogue/GAIA-DR1/Orig/';
MaxGet       = 10;


DataURL = sprintf('%s%s/',BaseURL,Type);

% Search for all URLS in the GAIA archive webpage
List    = www.find_urls(DataURL,'strfind',BaseFileName);

% retrieve all the GAIA files
PWD = pwd;
cd(CopyTo)
www.pwget(List,'',MaxGet)
cd(PWD)
