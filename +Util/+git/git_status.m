function [NewFiles,StrM,StrN]=git_status(varargin)
% Execute git status and return list of untrcaked files
% Package: Util
% Description: Execute git status and return list of untrcaked files.
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Show'      - Present status. Default is true.
%            'BasePath'  - Default is '/matlab/fun'.
% Output : - Cell array of relevant untracked files.
%          - Output string of modified files.
%          - Output string if untracked files.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [NewFiles,StrM,StrN]=Util.git.git_status
% Reliable: 2
%--------------------------------------------------------------------------


DefV.Show                 = true;
DefV.BasePath             = '/matlab/fun';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

InPar.BasePath = regexprep(InPar.BasePath,'/',filesep);

UserHome = Util.OS.get_userhome;

Path = sprintf('%s%s%s',UserHome,InPar.BasePath,filesep);
PWD = pwd;
cd(Path);


if (InPar.Show)
    [R,FullStr]=system('git status');
    fprintf('%s\n',FullStr)
end

[R,Str] = system('git -c color.status=false status --short --branch');

Tmp = regexp(Str,'\s+','split');
FlagQM = strcmp(Tmp,'??');
NewFiles = Tmp(find(FlagQM)+1);



% I1=strfind(Str,'modified content in submodules)')+31;
% I2=strfind(Str,'Untracked files:')-1;
% 
% StrM = Str(I1:I2);
% 
% I3=strfind(Str,'to include in what will be committed)')+37;
% I4=strfind(Str,'no changes added to commit (use "git add"')-1;
% 
% StrN = Str(I3:I4);
% Tmp = regexp(StrN,'\s','split');
% NewFiles = Tmp(~Util.cell.isempty_cell(Tmp));
Nn = numel(NewFiles);
Flag = false(Nn,1);
for In=1:1:Nn
    if (strcmp(NewFiles{In},'statrtup.m'))
        Flag(In) = true;
    else
        switch NewFiles{In}(1)
            case {'+','@'}
                Flag(In) = true;
            otherwise
                % do nothing
        end
    end
end
NewFiles = NewFiles(Flag);
            



cd(PWD);
