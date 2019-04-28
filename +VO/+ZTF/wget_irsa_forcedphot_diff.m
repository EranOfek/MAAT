function [Stat,Res]=wget_irsa_forcedphot_diff(RA,Dec,varargin)
% SHORT DESCRIPTION HERE
% Package: VO.ZTF
% Description: http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'User'    - String containing the IRSA/IPAC user name, or a
%                        a cell array containing a file name of user/pass.
%                        Default is
%                        {'/home/eran/matlab/passwords/ztf_ipac_pass'}.
%            'Pass'    - String containing the IRSA/IPAC password,.
%                        Default is [].
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

% RA = 280.8058788;
% Dec = 45.2077645;
% 
% RA = 241.518638;
% Dec = 36.871243;

OutFile = 'log.txt';

DefV.JDstart              = celestial.time.julday([1 9 2017]);
DefV.JDend                = celestial.time.julday([1 1 2021]);
DefV.User                 = {'/home/eran/matlab/passwords/ztfForced_ipac_pass'}; 
DefV.Pass                 = [];
DefV.email                = 'eran.ofek@weizmann.ac.il'; % note that in this service the e-mail is copuled to user/pass!
DefV.BaseURL              = 'http://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% get user/pass from passwords file
if (iscell(InPar.User))
    [InPar.User,InPar.Pass]=Util.files.read_user_pass_file(InPar.User{1});
end


URL = sprintf('%sra=%-10.6f&dec=%-10.6f&jdstart=%-11.3f&jdend=%-11.3f&email=%s',...
    InPar.BaseURL,RA,Dec,InPar.JDstart,InPar.JDend,InPar.email);
URL = replace(URL,' ','');

Options = weboptions('UserName',InPar.User,'Password',InPar.Pass);


webread(URL,Options);

CL = sprintf('wget --http-user=%s --http-passwd=%s -O %s "%s"',InPar.User,InPar.Pass,OutFile,URL);
[Stat,Res] = system(CL);
