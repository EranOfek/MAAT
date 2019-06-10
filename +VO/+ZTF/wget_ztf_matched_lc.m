function [DataTable,ColCell,Description]=wget_ztf_matched_lc(RA,Dec,varargin)
% wget ztf matched light curve data for sources around coordinates
% Package: VO.ZTF
% Description: wget ztf matched light curve data for sources around
%              coordinates.
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Band' - Default is 'g'.
%            'SearchRadius' - Default is 5.
%            'RadiusUnits'  = Default is 'arcsec'.
% Output : - Table of data.
%          - Cell array of column names.
%          - Cell array of column descriptions.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Data]=VO.ZTF.wget_ztf_matched_lc(298.0025, 29.87147)
% Reliable: 2
%--------------------------------------------------------------------------


DefV.Band                 = 'g';
DefV.SearchRadius         = 5;
DefV.RadiusUnits          = 'arcsec';
DefV.User                 = {'/home/eran/matlab/passwords/ztf_ipac_pass'}; 
DefV.Pass                 = [];
DefV.wgetProg             = 'wget';   % 'wget' | 'curl'
DefV.pwgetExtra           = '--no-check-certificate --load-cookies=cookies.txt';
DefV.CookiesFile          = 'cookies.txt';
DefV.OutTable             = 'ztf_lc_out.tbl';
DefV.AccountURL           = '/account/signon/login.do';
DefV.BaseURL              = 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

InPar.pwgetExtra = sprintf('%s -O %s',InPar.pwgetExtra,InPar.OutTable);

% get user/pass from passwords file
if (iscell(InPar.User))
    [InPar.User,InPar.Pass]=Util.files.read_user_pass_file(InPar.User{1});
end

% set up the IRSA cookies
[Stat,Res]=VO.ZTF.irsa_set_cookies('CookiesFile',InPar.CookiesFile,...
                                   'BaseURL',InPar.BaseURL,...
                                   'AccountURL',InPar.AccountURL,...
                                   'User',InPar.User,...
                                   'Pass',InPar.Pass,...
                                   'wgetProg',InPar.wgetProg);


SearchRadius = convert.angular(InPar.RadiusUnits,'deg',InPar.SearchRadius);                               
URL = sprintf('%s?POS=CIRCLE %f %f %f&BANDNAME=%s',InPar.BaseURL,RA,Dec,SearchRadius,InPar.Band);
% Example: https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE 298.0025 29.87147 0.0014&BANDNAME=g
%Data = webread(URL);
www.pwget(URL,InPar.pwgetExtra,1,[],'wget');

% read data from table
Data = Util.IO.file2str(InPar.OutTable,'str');

[DataTable,ColCell] = www.parse_html_table_data(Data);
Description = ColCell;
ColCell = {'ObjectID','ExposureID','HJD','MJD','Mag','MagErr','Flags',...
    'FilterCode','RA','Dec','Chi2','Sharpness','ExpFileTimeStamp',...
    'FieldID','CCDnumber','QuadrantID','LimMag5','MagZP','rmsMagZP',...
    'ColorTerm','ErrColorTerm','ExpTime','AirMass','ProgID'};
DataTable.Properties.VariableNames=ColCell;




                               