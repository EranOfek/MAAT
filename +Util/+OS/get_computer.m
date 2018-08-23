function Ans=get_computer
% Get computer name
% Package: Uti.OS
% Description: Get computer name.
% Input  : null
% Output : - Computer name.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Ans=Util.OS.get_computer
% Reliable: 2
%--------------------------------------------------------------------------

Ans = char(java.net.InetAddress.getLocalHost.getHostName);