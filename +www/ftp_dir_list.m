function File=ftp_dir_list(FTP,Recursive)
% Return files URLs from FTP containing a file listing.
% Package: www
% Description: Given an FTP URL that contains only a file listing,
%              return a cell array of all the files URLs.
% Input  : - FTP url.
%          - A flag {true|false} indicating if to recursively search for
%            files in all sub sirectories.
%            If false then will return only files and directories
%            in the current directory.
%            Default is true.
% Output : - Structure array of full URLs of files listed in FTP.
%            The structure contains the following fields:
%            .URL   - full URL
%            .isdir - is directory flag.
%            .name  - file name
%            .subdir- subdirectory relative to base
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FullURL=www.ftp_dir_list('ftp://legacy.gsfc.nasa.gov/chandra/data/science/ao01/cat1/1/primary/')
% Reliable: 2
%--------------------------------------------------------------------------
import www.*

if (nargin==1)
    Recursive = true;
end

Page = urlread(FTP);
if (isempty(Page))
    File = struct('URL',{},'isdir',{},'name',{},'subdir',{});
else
    SplitL  = regexp(Page,'\n','split');
    SplitS  = regexp(SplitL,'\s+','split');
    Nfile   = numel(SplitS);
    Ind = 0;
    for Ifile=1:1:Nfile
        if (numel(SplitS{Ifile})>1)
            Ind = Ind + 1;
            FileName = SplitS{Ifile}{end};
            File(Ind).URL    = sprintf('%s%s',FTP,FileName);
            File(Ind).name   = FileName;
            File(Ind).subdir = '';
            if (strcmp(SplitL{Ifile}(1),'d'))
                File(Ind).isdir = true;
                File(Ind).URL   = sprintf('%s/',File(Ind).URL);

                % clean '//' from URLs
                File(Ind).URL = regexprep(File(Ind).URL,'//','/');
                
                if (Recursive)
                    File(Ind).URL = sprintf('%s%s',File(Ind).URL,'/');
                    FileSubDir = ftp_dir_list(File(Ind).URL);
                    Nsubdir    = numel(FileSubDir);
                    File(Ind:Ind+Nsubdir-1) = FileSubDir;
                    [File(Ind:Ind+Nsubdir-1).subdir] = deal(FileName);
                    Ind = Ind+Nsubdir-1;

                end
            else
                File(Ind).isdir = false;
            end
        end
    end
end

