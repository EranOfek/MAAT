function [Cat,ColCell]=irsa_query_ztf_images(RA,Dec,varargin)
% Query ZTF images from the IRSA/IPAC database
% Package: VO.ZTF
% Description: Query ZTF images from the IRSA/IPAC database. Either raw,
%              science, cal images can be queried by position and metadata.
% Input  : - J2000.0 R.A. [radians, [H M S], or sexagesimal string], or
%            a string containing object name (e.g., 'm31').
%            If empty, then the query is done without position (only the
%            WHERE clause).
%          - J2000.0 R.A. [radians, [sign D M S], or sexagesimal string].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ImType'  - ZTF image type to query: 'sci' | 'raw' | 'cal'.
%                        Default is 'sci'.
%            'OutType' - Output table type: 'mat' | 'astcat'.
%                        Default is 'astcat'.
%            'Size'    - Position search size (heiht, [width]) [deg].
%                        Default is 0.
%            'SizeUnits'- Position search size units.
%                        Default is 'deg'.
%            'User'    - String containing the IRSA/IPAC user name, or a
%                        a cell array containing a file name of user/pass.
%                        Default is
%                        {'/home/eran/matlab/passwords/ztf_ipac_pass'}.
%            'Pass'    - String containing the IRSA/IPAC password,.
%                        Default is [].
%            'NameServer'- Name server function for interpreting the object
%                        name. Default is @VO.name.server_ned.
%            'Where'   - String containing WHERE clause (e.g.,
%                        'field=600 and ccdid=2'). Default is ''.
%            'Intersect'- Options are: 'COVERS' | 'ENCLOSED' | 'CENTER' |
%                        'OVERLAPS'.
%                        Default is 'OVERLAPS'.
%                        See https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html
%            'wgetProg'- Options are: 'wget' | 'curl'.
%                        Default is 'wget'.
%            'BaseURL' - Query website base URL:
%                        Default is 'https://irsa.ipac.caltech.edu'.
%            'AccountURL'- Default is '/account/signon/login.do'.
%            'DataURL'   - Default is '/ibe/search/ztf/products/'.
%            'CookiesFile' - Cookies file name 'cookies.txt'.
%            'OutTable'  - Output table. Default is 'out.tbl'.
%            'DelOut'  - Delete output table. Default is true.
% Output : - A matrix or an AstCat object containing the table of queried
%            images.
%          - A cell array of column names.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: T=VO.ZTF.irsa_query_ztf_images(358./RAD,23./RAD); % query by pos
%          T=VO.ZTF.irsa_query_ztf_images('m31');            % query ny name
%          T=VO.ZTF.irsa_query_ztf_images([],[],'Where','field=600 and ccdid=2');  
%          T=VO.ZTF.irsa_query_ztf_images([],[],'Where','obsjd>2458058 and obsjd<2458060');  
%          T=VO.ZTF.irsa_query_ztf_images([],[],'Where','field=600 and ccdid=2','ImType','raw'); 
% Reliable: 2
%--------------------------------------------------------------------------

CatField     = AstCat.CatField;
ColCellField = AstCat.ColCellField;

if (nargin<2)
    Dec = [];
    if (nargin<1)
        RA = [];
    end
end

    

DefV.ImType               = 'sci';  % 'sci' | 'raw' | 'cal'
DefV.OutType              = 'astcat';
DefV.Size                 = 0;
DefV.SizeUnits            = 'deg';
DefV.User                 = {'/home/eran/matlab/passwords/ztf_ipac_pass'}; 
DefV.Pass                 = [];
DefV.NameServer           = @VO.name.server_ned;
DefV.Where                = '';  % e.g., field=600
DefV.Intersect            = 'OVERLAPS'; % 'COVERS' | 'ENCLOSED' | 'CENTER' | 'OVERLAPS'
DefV.wgetProg             = 'wget';   % 'wget' | 'curl'
DefV.BaseURL              = 'https://irsa.ipac.caltech.edu';
DefV.AccountURL           = '/account/signon/login.do';
DefV.DataURL              = '/ibe/search/ztf/products/';
DefV.CookiesFile          = 'cookies.txt';
DefV.OutTable             = 'out.tbl';  % '-'
DefV.DelOut               = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% get user/pass from passwords file
if (iscell(InPar.User))
    [InPar.User,InPar.Pass]=Util.files.read_user_pass_file(InPar.User{1});
end

% query instructions:
% https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html
% https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html

% set up the IRSA cookies
[Stat,Res]=VO.ZTF.irsa_set_cookies('CookiesFile',InPar.CookiesFile,...
                                   'BaseURL',InPar.BaseURL,...
                                   'AccountURL',InPar.AccountURL,...
                                   'User',InPar.User,...
                                   'Pass',InPar.Pass,...
                                   'wgetProg',InPar.wgetProg);


Str.Query = '';
% construct query parameters
if (isempty(Dec))
    % Dec is empty
    % either RA contains object name
    % or RA is empty (only where search)
    if (isempty(RA))
        % WHERE search
        Str.POS = '';
    else
        % get RA/Dec from object name
        [RA,Dec]=InPar.NameServer(RA,'r');
        Str.POS = sprintf('POS=%f,%f&',RA,Dec);
    end
else
       
    RA  = celestial.coo.convertdms(RA,'gH','d');
    Dec = celestial.coo.convertdms(Dec,'gD','d');
    Str.POS = sprintf('POS=%f,%f&',RA,Dec);
end
if (isempty(Str.POS))
    Str.Size      = '';
    Str.Intersect = '';
else
    
    InPar.Size = convert.angular(InPar.SizeUnits,'deg',InPar.Size);
    InPar.Size = InPar.Size(:).*ones(2,1);

    Str.Size = sprintf('SIZE=%f,%f&',InPar.Size);

    Str.Intersect = sprintf('INTERSECT=%s&',InPar.Intersect);
end

if (isempty(InPar.Where))
    Str.Where = '';
else
    Str.Where = sprintf('WHERE=%s&',InPar.Where);
end

% %WHERE seeing  airmass  moonillf   maglimit  imgtype  obsjd fid  field ccdid
% %if (isempty(InPar.airmass))
% %    sprintf(
% Str.Where = '';
% %Str.Where = urlencode('WHERE=field=600');
% Str.Where = '&WHERE=field=600 AND ccdid=1';

Str.CT = sprintf('CT=csv');

%Str.Query = sprintf('%s&%s&%s&%s',Str.POS,Str.Size,Str.Intersect,Str.CT);
%Str.Query = sprintf('%s&%s',Str.Where,Str.CT);

Str.Query = sprintf('%s%s%s%s%s',Str.POS,Str.Size,Str.Intersect,Str.Where,Str.CT);

% run query
switch lower(InPar.wgetProg)
    case 'wget'
        [Stat, Res] = system(sprintf('wget --no-check-certificate --load-cookies=%s -q -O %s "%s%s%s?%s"',...
                        InPar.CookiesFile, InPar.OutTable, InPar.BaseURL, InPar.DataURL, InPar.ImType, Str.Query));
                    
        
    case 'curl'
        [Stat, Res] = system(sprintf('curl -b %s -o %s "%s%s%s?%s"',...
                        InPar.CookiesFile, InPar.OutTable, InPar.BaseURL, InPar.DataURL, InPar.ImType, Str.Query));

    otherwise
        error('Unknown wgetProg option');
end
if (Stat<0)
    disp(Res)
    error('Query failed');
end
%Res


[Table,ColCell] = VO.Util.read_casjobs_table(InPar.OutTable,'OutType','table');

switch lower(InPar.OutType)
    case 'mat'
        Cat = Table;
    case 'astcat'
        Cat = AstCat;
        Cat.(CatField) = Table;
        Cat.(ColCellField) = ColCell;
        Cat = colcell2col(Cat);
    otherwise
        error('Unknown OutType option');
end

if (InPar.DelOut)
    delete(InPar.OutTable);
    delete(sprintf('%s_1',InPar.OutTable));
end


% ZTF files location:
% https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html