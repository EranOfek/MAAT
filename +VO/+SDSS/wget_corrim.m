function LinkSt=wget_corrim(RA,Dec,varargin)
% Get corrected SDSS image
% Package: VO.SDSS
% Description: Retrieve SDSS corrected images by Run, Rerun, CamCol, Field
%              or RA/Dec.
%              Retrieve FITS files.
% Input  : - A 4 columns matrix of SDSS images IDs:
%            [Run, Rerun, CamCol, Field].
%            Or J2000.0 R.A. [rad, [H M S], sexagesimal
%            string].
%          - J2000.0 Dec. [rad, [Sign D M S], sexagesimal
%            string]. If empty, then assume first input
%            argument is [Run, Rerun, CamCol, Field] ID.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Band' - Bands of images to retrieve. Either a string
%                     containing a single or multiple bands
%                     (e.g., 'g' | 'ugr'), or a cell array of strings
%                     (e.g., {'g','r','i'}).
%            'Save' - Save the images. If false then only get the image
%                     links. If true then save the images in the current
%                     directory. If a dir name then cd to the dir name and
%                     retrieve the files.
%            'MaxGet'- Maximum number of parallel wget (see www.pwget.m).
%                     Default is 15.
%            'Extra' - Extra parameters to pass to wget (see www.pwget.m).
%                     Default is '-nc'.
%            'coo2runPar' - Cell array of parameters to pass to
%                      SDSS.coo2run.
%            'Unzip' - Unzip the images. Default is true.
% Output : - A structure array of images and links.
%            Link(ImageInd,BandInd)...
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Aug 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [LinkSt]=VO.SDSS.wget_corrim([1000 301 1 27]);
% Reliable: 2
%--------------------------------------------------------------------------

BaseURL = VO.SDSS.image_server; 

DefV.Band                = 'ugriz';
DefV.Save                = true;
DefV.MaxGet              = 15;
DefV.Extra               = '-nc --no-check-certificate';
DefV.coo2runPar          = {};
DefV.Unzip               = true;
DefV.GetFirst            = false;
%DefV.DelZip              = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (nargin==1)
    Dec = [];
end
if (isempty(Dec))
    % Assume RA contains ID
    ID = RA;
else
    % get ID from RA/Dec
    ID = cell2mat(VO.SDSS.coo2run(RA,Dec,InPar.coo2runPar{:}));
end

Nid = size(ID,1);

Run     = ID(:,1);
Rerun   = ID(:,2);
Camcol  = ID(:,3);
Field   = ID(:,4);

% Split a char array into cell array
Band = InPar.Band;
if (ischar(Band))
    TmpBand = Band;
    Nb      = numel(Band);
    Band    = cell(1,Nb); 
    for Ib=1:1:Nb
        Band{Ib} = TmpBand(Ib:Ib);
    end
end
Nband = numel(Band);        

%URL      = cell(Nid,1);
%FileName = cell(Nid,Nband);
%FullURL  = cell(Nid,Nband);

LinkSt   = Util.struct.struct_def({'FileName','Link','Band'},Nid,Nband);
for Iid=1:1:Nid
    % for each ID in ID matrix:
    URL      = sprintf('%s/%d/%d/%d/',BaseURL,Rerun(Iid),Run(Iid),Camcol(Iid));
    for Iband=1:1:Nband
        LinkSt(Iid,Iband).FileName = sprintf('frame-%s-%06d-%d-%04d.fits.bz2',Band{Iband},Run(Iid),Camcol(Iid),Field(Iid));
        LinkSt(Iid,Iband).Link     = sprintf('%s%s',URL,LinkSt(Iid,Iband).FileName);
        LinkSt(Iid,Iband).Band     = Band{Iband};
    end
end

PWD = pwd;
if (ischar(InPar.Save))
    % assume InPar.Save contains dir in which to save files
    cd(InPar.Save);
    InPar.Save = true;
end

% Retrieve the files
if (islogical(InPar.Save))
    if (InPar.Save)
        % get files
        %              After exceuting pwget.m it is difficult to kill it. In
        %              order to stop the execuation while it is running you
        %              have to create a file name 'kill_pwget' in the directory
        %              in which pwget is running (e.g., "touch kill_pwget").

        if (InPar.GetFirst && ~isempty(LinkSt))
            www.pwget({LinkSt(1).Link},InPar.Extra,InPar.MaxGet);
        else
            www.pwget({LinkSt(:).Link},InPar.Extra,InPar.MaxGet);
        end

        if (InPar.Unzip)
            % Unzip the images
            Nim = numel(LinkSt);
            for Iim=1:1:Nim
                system(sprintf('bunzip2 %s',LinkSt(Iim).FileName));
                FN = LinkSt(Iim).FileName(1:end-4);
%                             if (InPar.DelZip),
%                                 delete(LinkSt(Iim).FileName);
%                             end
                LinkSt(Iim).FileName = FN;
            end
        end

    end
else
    error('Illegal Save option');
end
% return to original dir
cd(PWD);


