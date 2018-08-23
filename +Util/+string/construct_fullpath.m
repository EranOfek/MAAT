function ProgFullPath=construct_fullpath(ProgName,ProgPath,CurFunDirName,CurFunName)
%--------------------------------------------------------------------------
% construct_fullpath function                                      General
% Description: Construct a full path string to a program name, given
%              its name, path or relative path and base path.
% Input  : - Program name.
%          - Program path. If startwith '.' then will be used as a relative
%            path (relative to the base path).
%            Otherwise, base path is ignored.
%          - Base path.
%          - Function name that resides in the base path. If given,
%            then will use which_dir.m to return base path.
% Output : - Program full path and name.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ProgFullPath=construct_fullpath('sex-2.5','../bin/SExtractor/src/','/home/eran/matlab/fun/ImPhot');
%          % or:
%          ProgFullPath=construct_fullpath('sex-2.5','../bin/SExtractor/src/',[],'sextractor');
%          % will return: '/home/eran/matlab/fun/ImPhot/../bin/SExtractor/src/sex-2.5
%          % and if absolute path is required:
%          ProgFullPath=construct_fullpath('sex-2.5','/bin/SExtractor/src/');
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==4)
    CurFunDirName = Util.files.which_dir(CurFunName);
end


if (isempty(ProgName))
    ProgName = '';
end

    % program path
    if (isempty(ProgPath))
        % assume program is in search path
        ProgFullPath = ProgName;
    else
        if (strcmp(ProgPath(end),filesep))
            FileSep = '';
        else
            FileSep = filesep;
        end

        if (strcmp(ProgPath(1),'.'))
            % use relative path
            ProgFullPath = sprintf('%s%s%s%s%s',CurFunDirName,filesep,ProgPath,FileSep,ProgName);
        else
            ProgFullPath = sprintf('%s%s%s',ProgPath,FileSep,ProgName);
        end
    end 

    
