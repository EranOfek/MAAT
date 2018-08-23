function Ans=get_user
% Get user name
% Package: Util.OS
% Description: Get user name.
% Input  : null.
% Output : - User name.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Ans=Util.OS.get_user;
% Reliable: 2
%--------------------------------------------------------------------------

Ans = char(java.lang.System.getProperty('user.name'));
