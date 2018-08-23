function [Sim,SimBO,OverScan,Data]=read_raw_image(List,varargin)
% Read ZTF raw image with its 4 quadrant and 4 bias overscan
% Package: ImUtil.ZTF
% Description: Read ZTF raw image with its 4 quadrant and 4 bias overscan.
%              Also construct a bias overscan image and subtract it from
%              the corresponding image.
% Input  : - A ZTF raw FITS image name of a list of images.
%            See Util.files.create_list for options.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'KeyQUADRANT          = 'QUADRANT';
%            'IndImages            = [2 3 4 5];
%            'IndOverScanImages    = [6 7 8 9];
%            'BiasOverscanBuffer   = [5 5];
%            'SigmaClip            = 4.5;
%            'MeanMethod           = @nanmean;
%            'MethodBias           = 'medfilt';
%            'MedFiltBlock         = 10;
%            'PolyDeg              = 2;
%            'GainCorrect          = true;
% Output : - A SIM object containing the bias overscan subtracted images.
%            The SIM object size is 4xN, where N is the number of images in
%            the input list and 4 is the number of amplifiers.
%          - A SIM object containing the corresponding bias overscan
%            images.
%          - A structure array of the same size as the first output
%            arguments containing the following fields:
%            'OverScan' - A vector of overscan values per row.
%            'Model'    - A smooth vector of overscan model values per row.
%          - A structure array containing the information breakout from the
%            image names.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Sim,SimBO]=ImUtil.ZTF.read_raw_image('ztf_20171107553194_000000_bi_c16_b.fits.fz');
%          [Sim,SimBO]=ImUtil.ZTF.read_raw_image('ztf_*_000000_bi_c16_b.fits.fz');
% Reliable: 2
%--------------------------------------------------------------------------

ImageField    = SIM.ImageField;
BackField     = SIM.BackField;
HeadField     = 'Header';


DefV.KeyQUADRANT          = 'QUADRANT';
DefV.IndImages            = [2 3 4 5];
DefV.IndOverScanImages    = [6 7 8 9];
DefV.BiasOverscanBuffer   = [5 5];
DefV.SigmaClip            = 4.5;
DefV.MeanMethod           = @nanmean;
DefV.MethodBias           = 'medfilt';
DefV.MedFiltBlock         = 10;
DefV.PolyDeg              = 2;
DefV.GainCorrect          = true;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% create cell arreay of list of images
[~,ListCell]=Util.files.create_list(List,NaN);

% number of images
Nim = numel(ListCell);
% number of parts (e.g., quadrants/amplifiers in each image)
Nparts = numel(InPar.IndImages);

% allocate memory for Sim images
% 4 x Nim
% 4 quadrant per image
Sim   = SIM(4,Nim);
% allocate memory for overscan Sim images
SimBO = SIM(4,Nim);

% read each image seperatly
OverScan = Util.struct.struct_def({'OverScan','Model'},Nparts,Nim);
for Iim=1:1:Nim
    % Read global array header
    [HeadGlobal,Nhdu]=FITS.get_head(ListCell{Iim},1);
    if (Nhdu~=9)
        error('Number of HDU in ZTF raw FITS image is not 9');
    end
    
    for Iparts=1:1:Nparts
        % read FITS image directly into SIM
        Sim(Iparts,Iim).(ImageField)   = FITS.read(ListCell{Iim},InPar.IndImages(Iparts));
        
        % read FITS overscan bias image directly into SIM
        SimBO(Iparts,Iim).(ImageField) = FITS.read(ListCell{Iim},InPar.IndOverScanImages(Iparts));
        
        % Read specific image header
        [HeadIm]=FITS.get_head(ListCell{Iim},InPar.IndImages(Iparts));
        % remove the block of keywords related to the fpack compression
        HeadIm.Header = HeadIm.Header(31:end,:);
        % combine global header and specific header
        Size = size(Sim(Iparts,Iim).(ImageField));
        HeadOp = add_key(HEAD,...
                              'NAXIS',2,'Number of dimensions',...
                              'NAXIS1',uint16(Size(2)),'Number of pixels in X-axis',...
                              'NAXIS2',uint16(Size(1)),'Number of pixels in Y-axis');
                           
        
        Head.(HeadField) = [HeadOp.(HeadField); HeadGlobal.(HeadField); HeadIm.(HeadField)];
                
        % Read header
        Sim(Iparts,Iim).(HeadField) = Head.(HeadField);      
        % Add amplifier keyword
        if (~isempty(InPar.KeyQUADRANT))
            Sim(Iparts,Iim) = add_key(Sim(Iparts,Iim),InPar.KeyQUADRANT,Iparts,'CCD quadrant/amplifier index');
        end
        
        %----------------------
        % overscan subtraction
        %----------------------
        BiasOverScan = SimBO(Iparts,Iim).(ImageField)(:,InPar.BiasOverscanBuffer(1):end-InPar.BiasOverscanBuffer(1));
        FlagPixUse   = double(abs(BiasOverScan - median(BiasOverScan,2)) < InPar.SigmaClip.*std(BiasOverScan,[],2));
        FlagPixUse(FlagPixUse==0) = NaN;
        MeanOverScan = InPar.MeanMethod(BiasOverScan.*FlagPixUse,2);
        X = (1:1:numel(MeanOverScan))';
        
        switch lower(InPar.MethodBias)
            case 'none'
                % do nothing
            case 'medfilt'
                % without the truncate option the median is biased at the
                % edges...
                
                BO = medfilt1(MeanOverScan,InPar.MedFiltBlock,'truncate');
                
                
            case 'polyfit'
                Par = polyfit(X,MeanOverScan,InPar.PolyDeg);
                BO  = polyval(Par,X);
                
            case 'row'
                BO = MeanOverScan;
                
            otherwise
                error('Unknown MethodBias option');
        end
        Sim(Iparts,Iim).(ImageField)  = Sim(Iparts,Iim).(ImageField) - BO;  
        % store bias overscan
        OverScan(Iparts,Iim).OverScan = MeanOverScan;
        OverScan(Iparts,Iim).Model    = BO;
         
        if (InPar.GainCorrect)
            Sim(Iparts,Iim) = gain_correct(Sim(Iparts,Iim));
        end
        
    end
    
    % optional stiching of pixels from other quadrants
    
    
    
end

if (nargout>2)
    Data = VO.ZTF.ztf_imagename_prop(ListCell);
end
