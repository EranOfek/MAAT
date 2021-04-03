function Created=cdmkdir(Dir)
% cd to directory - if not exist than create
% Package: +lastpipe.util
% Input  : - Directory name.
% Output : - A logical flag indicating if the directory was created (true)
%            or already exist (false).
% Example: lastpipe.util.cdmkdir('a/b')


try
    cd(Dir);
    Created = false;
catch
    if strcmp(Dir(1),filesep)
        cd(filesep);
    end
    Split = regexp(Dir,filesep,'split');
    % remove empty strings
    Split = Split(~Util.cell.isempty_cell(Split));
    Nsplit = numel(Split);
    for Isplit=1:1:Nsplit
        if exist(Split{Isplit},'dir')==0
            mkdir(Split{Isplit});
        end
        cd(Split{Isplit});
    end
    Created = true;
end


