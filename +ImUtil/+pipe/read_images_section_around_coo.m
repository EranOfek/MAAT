function [Sim,SimGlobalBack,Files,Section,CenterSection,IndGood]=read_images_section_around_coo(Files,RA,Dec,varargin)
% Read image sections around coordinates and select best images.
% Package: ImUtil
% Description: Read image sections around coordinates and select best
%              images. Images are read into SIM objects and by default 
%              the background and catalog fields are populated.
% Input  : - Cell array of file names, or file names with wild cards.
%            See Util.files.create_list.m for more options.
%          - J2000.0 R.A.
%          - J2000.0 Dec.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MexPar'  - Cell array of additional key,val parameters to
%                        pass to mextractor.m. Default is {}.
%            'ZP'      - mextractor ZP. Default is 22.
%            'MinNstar'- Minimum number of stars in image. Default is 20.
%                        If number of stars is smaller than remove image
%                        from list. If empty, then do not execute
%                        mextractor.
%            'MaxBackLevel'- Maximum background level.
%                        Default is 5000 (electrons).
%                        If background level is higher than remove image
%                        from list. If empty, then do not execute
%                        background.m.
%            'BackPar' - Cell array of additional key,val parameters to
%                        pass to background.m. Default is {}.
%            'SectionHalfSize'- Half size of image section to cut around
%                        coordinate. Default is [511 511].
%            'ExtraHalfSize'- Default is 3.
%            'MinDist' - Minimum distance of coordinate from image edge.
%                        Default is 64.
%            'Verbose' - Verbose. Default is true.
% Output : - Sim object containing all image sections.
%          - A SIM object with the global background values in the image
%            sections.
%          - Surviving file names.
%          - Section from original images [Xmin, Xmax, Ymin, Ymax].
%          - [X, Y] of center section in original images.
%          - Indices of good images returned in output.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Sim,SimGlobalBack,Files,Section,CenterSection]=ImUtil.pipe.read_images_section_around_coo(Files,RA,Dec)
% Reliable: 
%--------------------------------------------------------------------------

DefV.MexPar               = {};
DefV.ZP                   = 22;
DefV.MinNstar             = 20;
DefV.MaxBackLevel         = 5000;
DefV.BackPar              = {};
DefV.SectionHalfSize      = [511 511];
DefV.ExtraHalfSize        = 3;
DefV.MinDist              = 64;
DefV.Verbose              = true;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(Files))
    [~,Files] = Util.files.create_list(Files,NaN);
end

% Read header
Nim = numel(Files);
H   = FITS.get_head(Files);
W   = ClassWCS.populate(H);
% Get source X/Y position
X   = zeros(Nim,1);
Y   = zeros(Nim,1);
for Iim=1:1:Nim
    
    [X(Iim),Y(Iim)] = coo2xy(W(Iim),[RA,Dec]);
end
RoundXY = round([X,Y]);

% select overlapping X/Y regions that covers RA/Dec position.
ImageSizeXY = naxis(H);


% remove images with source near boundry
FlagInB = RoundXY(:,1) > InPar.MinDist & RoundXY(:,2)>InPar.MinDist & ...
         RoundXY(:,1) < (ImageSizeXY(1) - InPar.MinDist) & RoundXY(:,2) < (ImageSizeXY(2) - InPar.MinDist);

Files = Files(FlagInB);
ImageSizeXY = ImageSizeXY(FlagInB,:);
RoundXY = RoundXY(FlagInB,:);

%[Section,FlagIn,CenterSection] = ImUtil.Im.best_section(ImageSizeXY, RoundXY, InPar.SectionHalfSize+InPar.ExtraHalfSize, InPar.MinDist);
[Section,FlagIn,CenterSection] = ImUtil.Im.adaptive_section(ImageSizeXY, RoundXY, InPar.SectionHalfSize+InPar.ExtraHalfSize, InPar.MinDist);

% remove images in which the requested RA/Dec is near boundry
Files = Files(FlagIn);
Section = Section(FlagIn,:);
Nim = numel(Files);

FlagInB(find(FlagInB)) = FlagIn;


if (InPar.Verbose)
    fprintf('%d images were selected out of %d\n',Nim,numel(FlagIn));
end

% Read image sections to SIM object
Sim = FITS.read2sim(Files,'CCDSEC',Section);

% correct for gain
Sim = gain_correct(Sim);

SimGlobalBack = [];
IndGood       = true(size(Sim));
if ~isempty(InPar.MaxBackLevel)
    % Sim local background
    Sim = background(Sim,InPar.BackPar{:});
    % Sim globacl background
    SimGlobalBack = background(Sim,'Block','full');

    % remove images with very high background level
    FlagGoodBack  = [SimGlobalBack.BackIm]<InPar.MaxBackLevel;
    IndGood       = find(FlagGoodBack);
    % select only images with good background level
    SimGlobalBack = SimGlobalBack(FlagGoodBack);
    Sim           = Sim(FlagGoodBack);
    Files         = Files(FlagGoodBack);
    Section       = Section(FlagGoodBack,:);
end

if ~isempty(InPar.MinNstar)
    % Extract sources for alignment ---
    [Sim,MetaData] = mextractor(Sim,'ZP',InPar.ZP,InPar.MexPar{:});

    Ncat = sizecat(Sim);
    FlagCat = Ncat>InPar.MinNstar;
    IndGood = IndGood(FlagCat);
    SimGlobalBack = SimGlobalBack(FlagCat);
    Sim           = Sim(FlagCat);
    Files         = Files(FlagCat);
    Section       = Section(FlagCat,:);
end
        