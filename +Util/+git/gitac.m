function [Res,Stat]=gitac(List,Comment)
% Add a list of files to a git repository and commit.
% Package: www
% Description: add a list of files to a git repository and commit.
% Input  : - A file name or a cell array of file names to add and commit.
%            The file path will be determined using the which command.
%          - An optional comment. Default is ''
% Output : - Result.
%          - Status.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: www.gitac('gitac.m');
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
    Comment = '';
end

if (iscell(List))
    Ncell = numel(List);
    for Icell=1:1:Ncell
        [Stat,Res]=Util.git.git('add',which(List{Icell}));
    end
elseif (ischar(List))
    [Stat,Res]=Util.git.git('add',which(List));
else
    error('Invalid List type');
end

% commit
Comment = sprintf('"%s"',Comment);
[Stat,Res] = Util.git.git('commit','-m',Comment);


