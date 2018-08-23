function [Stat,Res]=irsa_set_cookies(varargin)
% Set user/pass cookies for IRSA query
% Package: VO.ZTF
% Description: Set user/pass cookies for IRSA query
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'User'    - String containing the IRSA/IPAC user name, or a
%                        a cell array containing a file name of user/pass.
%                        Default is
%                        {'/home/eran/matlab/passwords/ztf_ipac_pass'}.
%            'Pass'    - String containing the IRSA/IPAC password,.
%                        Default is [].
%            'wgetProg'- Options are: 'wget' | 'curl'.
%                        Default is 'wget'.
%            'BaseURL' - Query website base URL:
%                        Default is 'https://irsa.ipac.caltech.edu'.
%            'AccountURL'- Default is '/account/signon/login.do'.
%            'CookiesFile' - Cookies file name 'cookies.txt'.
% Output : - Status.
%          - Result string.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Stat,Res]=VO.ZTF.irsa_set_cookies;
% Reliable: 
%--------------------------------------------------------------------------

DefV.User                 = {'/home/eran/matlab/passwords/ztf_ipac_pass'}; 
DefV.Pass                 = [];
DefV.wgetProg             = 'wget';   % 'wget' | 'curl'
DefV.BaseURL              = 'https://irsa.ipac.caltech.edu';
DefV.AccountURL           = '/account/signon/login.do';
DefV.CookiesFile          = 'cookies.txt';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% get user/pass from passwords file
if (iscell(InPar.User))
    [InPar.User,InPar.Pass]=Util.files.read_user_pass_file(InPar.User{1});
end


% query instructions:
% https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html

% set up the IRSA cookies
switch lower(InPar.wgetProg)
    case 'wget'
        [Stat,Res] = system(sprintf('wget --no-check-certificate --save-cookies=%s -q -O /dev/null "%s%s?josso_cmd=login&josso_username=%s&josso_password=%s"',...
                        InPar.CookiesFile, InPar.BaseURL, InPar.AccountURL, InPar.User, InPar.Pass));
    case 'curl'
        [Stat, Res] = systen(sprintf('curl -c %s "%s%s?josso_cmd=login&josso_username=%s&josso_password=%s"',...
                        InPar.CookiesFile, InPar.BaseURL, InPar.AccountURL, InPar.User, InPar.Pass));
    otherwise
        error('Unknown wgetProg option');
end
if (Stat<0)
    disp(Res)
    error('Set cookies failed');
end
