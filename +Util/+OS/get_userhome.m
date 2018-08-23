function Ans=get_userhome
% Get user home directory path
% Package: Util.OS
% Description: Get user home directory path.
% Input  : null.
% Output : - User home directory path.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Ans=Util.OS.get_userhome
% Reliable: 2
%--------------------------------------------------------------------------

Ans   = char(java.lang.System.getProperty('user.home'));
