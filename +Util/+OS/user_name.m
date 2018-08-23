function User=user_name
% Get user name
% Package: Util.OS
% Description: Get the current user name.
% Input  : null
% Output : - User name.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%   BUGS : Does not work on Windows
% Example: User=Util.OS.user_name;
% Reliable: 2
%--------------------------------------------------------------------------


if (isunix)
    [~,User] = system('echo $USER');
    User     = Util.string.spacedel(User);
end

if (ismac)
    [~,User] = system('echo $USER');
    User     = Util.string.spacedel(User);
end

if (ispc)
    error('does not work on Windows');
end

