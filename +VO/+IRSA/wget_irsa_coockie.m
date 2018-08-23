function CookieName=wget_irsa_coockie(varargin)
% Get IRSA cookie for a user and password
% Package: VO.IRSA
% Description: Get IRSA cookie for a user and password
% Input  : - User name. Default is 
%          - Password. Default is
%          - Cookie name. Default is 'ptf.txt'.
% Output : - String containing cookie name.
% Example: CookieName=VO.IRSA.wget_irsa_coockie
% Reliable: 
%--------------------------------------------------------------------------

Def = {'PTF','palomar','ptf.txt'};
Ndef = numel(Def);

Narg = numel(varargin);
Arg  = [varargin, Def(Narg+1:Ndef)];
User = Arg{1};
Pass = Arg{2};
CookieName = Arg{3};


Command = sprintf('wget --save-cookies=%s -O /dev/null "http://irsa.ipac.caltech.edu/account/signon/login.do?josso_cmd=login&josso_username=%s&josso_password=%s"',CookieName,User,Pass);
system(Command);

