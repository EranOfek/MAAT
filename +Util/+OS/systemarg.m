function [Res,Stat]=systemarg(varargin)
% Running the UNIX system command.
% Package: Util.OS
% Description: Running the UNIX system command.
% Input  : * Arbitrary number of arguments that will be concatenated
%            with spaces to a single string and will be run using
%            the system command.
%            If you want to use ' sign use ''' instead.
% Output : - Results.
%          - Status.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: import Util.OS.*
%          systemarg sed "s#_# #g" try > try1
%          systemarg awk '''{ print }''' /etc/passwd
% Reliable: 2
%--------------------------------------------------------------------------

Narg = numel(varargin);
Cmd = '';
for Iarg=1:1:Narg
    Cmd = sprintf('%s %s',Cmd,varargin{Iarg});
end

[Stat,Res] = system(Cmd);
