function [User,Pass]=read_user_pass_file(UPfile)
% Read user/password from file
% Package: Util.code
% Description: Read user/password from file
% Input  : - File name and path.
% Output : - User
%          - Password
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [User,Pass]=Util.files.read_user_pass_file('~/matlab/passwords/ztf_archive_pass');
% Reliable: 2
%--------------------------------------------------------------------------


% assume user/pass are given in file
FID = fopen(UPfile,'r');
User = fgetl(FID);
Pass = fgetl(FID);
fclose(FID);
