function [DirsList,DirsListURL,FilesList,FilesListURL]=rftpget(URL,varargin)
% Recursively retrieve the entire directory tree in an FTP site
% Package: www
% Description: A wrapper around the wget command designed to recursively
%              retrieve the entire directory tree in an FTP site.
%              See also: mwget.m, pwget.m, find_urls.m
% Input  : - A string containing a URL.
%            Alternativly this can be a cell array of URLs.
%            However, if this is a cell array than 'rftp' will be set to 0.
%          * Arbitary number of ...,keyword,value,... pairs.
%            Where allowed keywords are:
%            'user' - user
%            'pass' - password
%            'par'  - string of additional parameters to pass to wget.
%            'out'  - Output file name
%            'rftp' - recursive wget of an FTP site (will get all URLs
%                     for all files under a FTP directory.
%                     keyword is 0 | 1, default is 0 (no).
%                     Note that this option will get only files names,
%                     but not the files themselfs.
%            'hidden' - how to deal with hidden files/directories in FTP
%                     structure. If 0, then don't retrieve hidden files,
%                     if 1 than retrieve hidden files. Default is 0.
%            'rmindex' - rm index.html files at the end if process.
%                     0 - don't remove; 1 - remove (default).
% Output : - In case the URL is an ftp directory return all the 
%            directories listed in this URL.
%          - In case the URL is an ftp directory return all the 
%            full URLs of the directories listed in this URL.
%          - In case the URL is an ftp directory return all the 
%            files listed in this URL.
%          - In case the URL is an ftp directory return all the 
%            full URLs of the files listed in this URL.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    May 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DirsList,DirsListURL,FilesList,FilesListURL]=www.rftpget('ftp://legacy.gsfc.nasa.gov/chandra/data/science/ao01/cat1/1/','rftp',1);
% Reliable: 2
%--------------------------------------------------------------------------
FTP_Listing = '.listing';
IsFTP = ~isempty(strfind(URL,'ftp://'));


FilesList    = {};
FilesListURL = {};
DirsList     = {};
DirsListURL  = {};

Par     = '';
RecFTP  = 0;
Hidden  = 0;
RmIndex = 1;
AddArg = ' ';   % string of additional arguments to pass to wget
Narg = length(varargin);
for Iarg=1:2:Narg-1,
   switch lower(varargin{Iarg})
    case 'user'
       if (IsFTP),
	  AddArg = sprintf('%s--ftp-user=%s ',AddArg,varargin{Iarg+1});
       else
	  AddArg = sprintf('%s--user=%s ',AddArg,varargin{Iarg+1});
       end
    case {'pass','password'}
       if (IsFTP),
	  AddArg = sprintf('%s--ftp-upassword=%s ',AddArg,varargin{Iarg+1});
       else
	  AddArg = sprintf('%s--upassword=%s ',AddArg,varargin{Iarg+1});
       end
    case 'par'
       Par = sprintf(' %s',varargin{Iarg+1});
    case {'o','out'}
       AddArg = sprintf('%s--output-document=%s ',AddArg,varargin{Iarg+1});
    case 'rftp'
       RecFTP = varargin{Iarg+1};
    case 'hidden'
       Hidden = varargin{Iarg+1};
    case 'rmindex'
       RmIndex = varargin{Iarg+1};
    otherwise
        error('Unknown keywords option');
   end
end

% read .listing file from ftp servers:
if (IsFTP),
   AddArg = sprintf('%s--no-remove-listing ',AddArg);
end



if (nargout>0),
   if (IsFTP),
      ListingFileCell = Util.files.file2str(FTP_Listing,'cell');
      FilesInDir      = regexp(ListingFileCell,' ','split');
   
      Nd = length(FilesInDir);
      Idir = 0;
      Ifile = 0;
      for Id=3:1:Nd,
         if (strcmp(FilesInDir{Id}{end}(1),'.') & Hidden==0),
            % dont retrieve hidden file
         else
            if strcmp(FilesInDir{Id}{1}(1),'d'),
               % directory
               Idir = Idir + 1;
               DirsList{Idir}     = FilesInDir{Id}{end};
               DirsListURL{Idir}  = sprintf('%s%s%s',URL,DirsList{Idir},'/');
            else
               % file
               Ifile = Ifile + 1;
               FilesList{Ifile}    = FilesInDir{Id}{end};
               FilesListURL{Ifile} = sprintf('%s%s',URL,FilesList{Ifile});
            end
         end
      end
   
      delete(FTP_Listing);
  
   
      switch RecFTP
       case 1
          %--- recursive FTP ---
          Ndir = length(DirsList);
          for Idir=1:1:Ndir,
             NewURL = sprintf('%s%s%s',URL,DirsList{Idir},'/');
             [NewDirsList,NewDirsListURL,NewFilesList,NewFilesListURL] = mwget(NewURL,varargin{:});
             DirsList     = [DirsList, NewDirsList];
             DirsListURL  = [DirsListURL, NewDirsListURL];
             FilesList    = [FilesList, NewFilesList];
             FilesListURL = [FilesListURL, NewFilesListURL];
          end
       otherwise
          % do nothing
      end
   end
end
   
% remove index.html files
switch RmIndex
 case 1
    delete('index.htm*');
 otherwise
    % do nothing
end
   
