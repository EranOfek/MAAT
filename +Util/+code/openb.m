function openb(FileName)
% Open matlab editor and save a backup copy of previous file
% Package: Util.code
% Description: open a matlab function in matlab editor (like the open
%              command). Before the file is opened it is backuped in
%              a backup (default is 'old') directory in the file directory
%              under the name FileName.current_date
% Input  : - File name to backup and open.
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

BackupDir = 'old';
Path = which_dir(FileName);

if (~strcmp(FileName(end-1:end),'.m'))
    FileName = sprintf('%s.m',FileName);
end

PWD = pwd;
cd(Path);
if ~isdir(BackupDir)
    mkdir(BackupDir);
end
Time = clock;
copyfile(FileName,sprintf('.%s%s%s%s.%02d%02d%04d',filesep,BackupDir,filesep,FileName,Time(3),Time(2),Time(1)));
cd(PWD);

open(FileName);
