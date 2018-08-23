function URL=wget_obsid(ObsID,varargin)
% Get all the files associated with a Chandra ObsID
% Package: VO.Chandra
% Description: Get all the files associated with a Chandra ObsID
%              The Chandra observations catalog is stored in
%              cats.X.ChandraObs.
%              Use VO.Chandra.build_obsid_cat to build catalog.
% Input  : - ObsID
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Download' - Indicating what to download:
%                         'all' - all data (default).
%                         'primary' - only data in primary directory.
%                         'secondary' - only data in secondary directory.
%            'Output'   - Output data in a 'flat' structure' or
%                         'dir' (directory) structure.
%                         'none' - do not copy.
%                         Default is 'dir'.
%            'Ungzip'   - ungzip data {true|false}. Default is true.
%            'CopyTo'   - Directory in which to cd before copying data.
%            'ReGet'    - Get files again if 'primary' dir exist 
%                         in directory {true|false}. Default is false.
%            'Extra'    - Extra parameters for wget.
%            'MaxGet'   - Maximum number of files for parallel get
%                         (see www.pwget.m). Default is 10.
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.Chandra.wget_obsid(366)
% Reliable: 2
%--------------------------------------------------------------------------
import Util.IO.*

%DefV.ChandraCat    = [];
DefV.Download      = 'all';   % {'all'|'primary'|'secondary'}
DefV.Output        = 'dir';   % {'dir'|'flat'}
DefV.Ungzip        = true;
DefV.CopyTo        = [];
DefV.ReGet         = false;
DefV.Extra         = '';
DefV.MaxGet        = 10;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% if (isempty(InPar.ChandraCat))
%     Cat = VO.Chandra.build_obsid_cat('GetInfo',false,'Verbose',false);
% else
%     if (ischar(InPar.ChandraCat))
%         Cat = load2(InPar.ChandraCat);
%     else
%         Cat = InPar.ChandraCat;
%     end
% end

Cat = cats.X.ChandraObs;


% Search ObsID
%Iobs = find([Cat.ObsID]==ObsID);
Iobs  = find(Cat.Cat.ObsID==ObsID);
%URL  = Cat(Iobs).URL;
URL   = Cat.Cat.URL{Iobs};



PWD = pwd;
if (InPar.CopyTo)
    cd(InPar.CopyTo);
end

% check if directory is populated
if (exist('primary','dir')>0 && ~InPar.ReGet)
    % do not get files again
else

    % get file names
    switch lower(InPar.Download)
        case 'all'
            List = www.ftp_dir_list(URL);
            I1 = find(strcmp({List.subdir},''));
            Ip = find(strcmp({List.subdir},'primary'));
            Is = find(strcmp({List.subdir},'secondary'));
        case 'primary'
            List = www.ftp_dir_list(sprintf('%sprimary/',URL));
            I1   = [];
            Ip   = (1:1:length(List))';
            Is   = [];
        case 'secondary'
            List = www.ftp_dir_list(sprintf('%ssecondary/',URL));
            I1   = [];
            Ip   = [];
            Is   = (1:1:length(List))';
        otherwise
            error('Unknown Download option');
    end


    % download
    switch lower(InPar.Output)
        case 'none'
            % do nothing
            
        case 'flat'
            www.pwget({List.URL},InPar.Extra,InPar.MaxGet);        
            if (InPar.Ungzip)
                system('gzip -d *.gz');
            end
        case 'dir'
            www.pwget({List(I1).URL},InPar.Extra,InPar.MaxGet);
            if (~isempty(Ip))
                mkdir('primary');
                cd('primary');
                www.pwget({List(Ip).URL},InPar.Extra,InPar.MaxGet);
                if (InPar.Ungzip)
                    system('gzip -d *.gz');
                end
                cd ..
            end
            if (~isempty(Is))
                mkdir('secondary');
                cd('secondary');
                www.pwget({List(Is).URL},InPar.Extra,InPar.MaxGet);
                if (InPar.Ungzip)
                    system('gzip -d *.gz');
                end
                cd ..
            end
       otherwise
            error('Unknown Output option');
    end

end  
cd(PWD);                        