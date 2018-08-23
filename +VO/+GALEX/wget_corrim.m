function [LinkSt]=wget_corrim(RA,Dec,varargin)
% Get GALEX corrected images from image server
% Package: VO.GALEX
% Description: Get GALEX corrected images
% Input  : - J2000.0 R.A. [rad, [H M S], or sexagesimal
%            string], or and Index in the GALEX images file.
%          - J2000.0 Dec. [rad, [Sign D M S], or sexagesimal
%            string]. If first argument is ID then this need to
%            be an empty matrix. Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Band' - Bands of images to retrieve. Either a string
%                     containing a single or multiple bands
%                     (e.g., 'n' | 'nf'), or a cell array of strings
%                     (e.g., {'n,'f'}).
%            'Save' - Save the images. If false then only get the image
%                     links. If true then save the images in the current
%                     directory. If a dir name then cd to the dir name and
%                     retrieve the files.
%            'ImType'- A string or cell array of image types to retrieve.
%                      'i'-int.
%                      'c'-cnt.
%                      'r'-rrhr(exp time).
%                      Default is 'icr'.
%            'DR'    - Data release. Default is 'GR6'.
%            'MaxGet'- Maximum number of parallel wget (see www.pwget.m).
%                     Default is 15.
%            'Extra' - Extra parameters to pass to wget (see www.pwget.m).
%                     Default is '-nc'.
%            'Unzip' - Unzip the images. Default is true.
% Output : - Structure array of links to GALEX images.
%            Link(ImageInd,BandIndex,ImTypeIndex)...
% Example: [LinkSt]=VO.GALEX.wget_corrim(0,0.1);
% Reliable: 2

if (nargin==1)
    Dec = [];
end

GALEX_Server = VO.GALEX.image_server;

DefV.Band                = 'fn';
DefV.Save                = true;
DefV.ImType              = 'icr';   % i-intensity, c-counts, r-rrhr(exp)
DefV.DR                  = 'GR7';
DefV.MaxGet              = 15;
DefV.Extra               = '-nc';
DefV.Unzip               = true;
%DefV.DelZip              = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (isempty(Dec))
    % assume RA contains images ID
    Ind = RA;
else
    [~,~,Ind] = VO.GALEX.coo2id(RA,Dec,InPar.DR);
end
% convert ID to matrix
if (iscell(Ind))
    Ind = cell2mat(Ind);
end
Ind = Ind(:);
[PathNUV,PathFUV]=VO.GALEX.ind2path(Ind,InPar.DR);

% Filter
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

% ImType
% Split a char array into cell array
ImType = InPar.ImType;
if (ischar(ImType))
    TmpImType = ImType;
    Nt      = numel(ImType);
    ImType  = cell(1,Nt); 
    for It=1:1:Nt
        ImType{It} = TmpImType(It:It);
    end
end
Ntype = numel(ImType);   

% number of image fields
Nim = numel(Ind);

% for each image field

LinkSt = Util.struct.struct_def({'FileName','Link','Band','ImageType'},0,0);

for Iim=1:1:Nim
    % for each image

    for Iband=1:1:Nband
        % for each band
        Filter = Band{Iband};
        switch lower(Filter)
            case 'f'
                FilePath = PathFUV{Iim};
            case 'n'
                FilePath = PathNUV{Iim};
            otherwise
                error('Unknown Filter option');
        end
        if (~any(isnan(FilePath)) && ~strcmp(FilePath,'NaN'))
            for Itype=1:1:Ntype
                % for each image type
                switch lower(ImType{Itype})
                    case 'i'
                        ImageType = 'int';
                    case 'c'
                        ImageType = 'cnt';
                    case 'r'
                        ImageType = 'rrhr';
                    otherwise
                        error('Unknown ImType option');
                end
                FilePathType = regexprep(FilePath,'-int.',sprintf('-%s.',ImageType));
                Splitted = regexp(FilePathType,'/','split');

                LinkSt(Iim,Iband,Itype).FileName  = Splitted{end};
                LinkSt(Iim,Iband,Itype).Link      = sprintf('%s%s',GALEX_Server,FilePathType);
                LinkSt(Iim,Iband,Itype).Band      = Filter;
                LinkSt(Iim,Iband,Itype).ImageType = ImageType;
            end
        end

    end
end

% Save
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
        %              After exceuting www.pwget.m it is difficult to kill it. In
        %              order to stop the execuation while it is running you
        %              have to create a file name 'kill_pwget' in the directory
        %              in which pwget is running (e.g., "touch kill_pwget").
        www.pwget({LinkSt.Link},InPar.Extra,InPar.MaxGet);

        if (InPar.Unzip)
            % Unzip the images
            Nim = numel(LinkSt);
            for Iim=1:1:Nim
%                             Iim
%                             LinkSt(Iim).FileName
                %FN = gunzip(LinkSt(Iim).FileName);
                system(sprintf('gunzip %s',LinkSt(Iim).FileName));
                %LinkSt(Iim).FileName
                %Status
                FN = LinkSt(Iim).FileName(1:end-3);
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




            