function [Res,Stat]=git_cmd(varargin)
% Call the git command with parameters
% Package: www
% Description: Call the git command with parameters
% Input  : * Arbitrary number of arguments to pass to the git command.
% Output : - Result.
%          - Status.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: gitac.m
% Example: Util.git.git_cmd add *.m
%          Util.git.www.git_cmd commit -m "initial project version"
% Reliable: 2
%--------------------------------------------------------------------------

Narg = numel(varargin);

Str = '';
for Iarg=1:1:Narg
    Str = sprintf('%s %s',Str,varargin{Iarg});
end

[Stat,Res] = system(sprintf('git %s',Str));
