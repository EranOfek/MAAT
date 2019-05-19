function [Cat,MetaData]=mextractor(Sim,varargin)
% Matched filter source extraction and measurments
% Package: @SIM
% Description: Source extraction, photometry and measurments using matched
%              filters.
%              Search for sources simultanously using several matched
%              filters (e.g., the PSF and exponential profiles), threshold
%              the matched filter image, search for peaks and measure the
%              properties of the peaks.
%              Main features:
%              - Automatic estimation of the PSF (2 pass mode). This run
%              the function twice. First to look for bright PSF stars, and
%              next to find and measure all the sources.
%              - Measuring the properties of sources (only or also) in a
%              pre-defined list of position (Force position parameters).
%              - Optimal search for point sources and optional search for
%              extended sources (simultanously).
%              - Measuring PSF photometry.
%              - Measuring aperture photometry.
%              - Measure background locally and globally.
%              - Measuring S/N in the filtered images and unfilter image.
%              - Calculate the airmass, parallactic angle, azimuth,
%              altitude for each source.
%              - Given a mask image, propagate the bit mask for each source
%              into the source FLAGS.
%              - Return the output magnitude in either magnitudes or
%              luptitudes.
% Input  : - A SIM object. The SIM object should contain an image in the
%            image field. Optionally it may contain a background image, and
%            std image, a mask image and a PSF. If these extra field are
%            not available they will be calculated and populated.
%            See images2sim.m for reading images into SIM.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Output properties ---
%            'ColCell'  - A cell array of column names to include in the
%                         output. The following column names are available:
%                         'NUMBER'
%                         'ISFORCE'
%                         
%            'ColFile'  - A file name 
%                         containing in the first column the requested
%                         output column names (listed in the 'ColCell'
%                         argument description.
%                         An example for such a file is given in
%                         'mextractor.outpar'.
%                         If empty then use the 'ColCell' argument.
%                         If given, this parameter override the 'ColCell'
%                         argument. Default is empty.
%            --- Primary detection parameters ---
%            'Thresh' - Detection threshold in units of sigma (filtered
%                       image std). Default is 5.
%            --- Secondary detection parameters (NO NEED TO MODIFY) --- 
%            'MinArea'     - Minimum area of source above threshold.
%                            Default is 1.
%            'AreaOpenConn'- Connectivity paramaeter for threshold.m.
%                            See bwareaopen.m for description.
%                            Default is 8.
%            'RegionMaxConn'- Connectivity parameter for imregionalmax.m
%                            (used in threshold.m and local_maxima.m).
%                            Default is 8.
%            'ReplaceVal'   - If not empty then replace NaN values in the
%                            MaxImage by this scalar.
%                            See threshold.m for details.
%                            Default is [].
%            'ColXY'        - Output position columns for threshold.m.
%                            Default is {'XPIX_PEAK','YPIX_PEAK'}.
%            --- Forced position ---
%            'ForcePos'     - Forced positions in which to measure source
%                            properties, even if a source was not found.
%                            This is a two column matrix of [X,Y] position
%                            or an AstCat object containing the [X,Y]
%                            columns.
%                            Default is [] (i.e., no force position).
%            'ForceCatCol'  - If 'ForcePos' is an AstCat object then this
%                            is a cell array of column names containing
%                            the X,Y positions at which to force measure
%                            the sources.
%                            Default is {'XWIN_IMAGE','YWIN_IMAGE'};
%            'OnlyForce'    - A logical flag indicating if to measure only
%                            the force position list (true) or to
%                            search for sources and measure both the found
%                            sources and the forced list (false).
%                            Default is false.
%            'WinPosFromForce' - A logical flag indicating if the output
%                            windowed positions are the one provided in
%                            the forced position list (true), or should
%                            be re-calculated by the function (false).
%                            Default is true.
%            --- Filters ---
%            'FilterFind'  - A logical flag (true|false) indicating if
%                            to run mextractor in double-pass mode. In this
%                            mode mextractor first find bright stars using
%                            whatever default filter the user provided and
%                            use them to construct the image PSF. In the
%                            second pass mextractor executed with the PSF
%                            found in the 1st iteration as the filter.
%                            Default is true.
%                            Note that if the 'Filter' argument is provided 
%                            or the PSF is included in the SIM then this
%                            will be effectively set to false.
%            'FilterFindThresh' - Threshold, in units of sigmas, for
%                            finding bright stars for PSF extraction.
%                            Default is 10.
%            'FilterIm'    - A logical flag indicating if to filter the
%                            image prior to thresholding. Default is true.
%            'Filter'             = []; %@Kernel2.gauss;
%            'FilterFunPar'       = {1.5,1.5,0,15,15};
%            'DefFilter'          = @Kernel2.gauss;
%            'DefFilterFunPar'    = {1.5,1.5,0,15,15};
%            'FilterFun'          = 'auto';
%            'GetPsfPar'          = {};    % additional parameters to pass to getmpsf.m
%            'AddFilter'          = {@Kernel2.exp};    % cell array of filters to use on all images
%            'AddFilterPar'       = {{3,15,15}};
%            'ConvAddFilter'      = true;
%            --- Magnitudes and zero points ---
%            'ZP'                 = 22;
%            'MagInLuptitude'     = true;
%            'LuptSoft'           = 1e-10;
%            --- PSF photometry ---
%            'PsfFitRad'          = 2;   % pix
%            --- Aperture photometry ---
%            'AperRad'            = [2 4 8 12 16];
%            'Annulus'            = [20 24];
%            'AperBackFun'        = @median; %@rmean; %@median;
%            'AperBackFunPar'     = {}; %{1,[0.25 0.25]};  %{};
%            'AperBackErrFun'     = @std;
%            'AperBackErrFunPar'  = {};
%            --- Gain ---
%            'Gain'               = {'GAIN'};  % or numeric
%            'OrigGainKey'        = 'ORIGGAIN';
%            --- Background ---
%            'BackSub'            = false;   % indicate if back already subtracted
%            'ReBackStd'          = false;
%            'BackPar'            = {};
%            '
% Output : - The input SIM object were the image is in units of electrobs
%            (i.e., multiplied by the gain), or a new AstCat object with the
%            extracted source catalog populated in the 'Cat' field.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: SL=[500 500 1e5; 200.12 203.89 1e5];
%          Image=image_art('ImSize',[1024 1024],'StarList',SL,'Back',300)
% Reliable: 2
%--------------------------------------------------------------------------
SN2DM = 2.5./log(10);  % S/N to mag error factor - 1.086...
RAD   = 180./pi;

ImageField     = SIM.ImageField;
BackField      = SIM.BackField;
ErrField       = SIM.ErrField;
CatField       = 'Cat';
ColCellField   = 'ColCell';
MaskField      = 'Mask';

%--------------------------
%--- Default parameters ---
%--------------------------
% Detection parameters
DefV.Thresh             = 5;      % sigmas
DefV.ThreshIsSigma      = true;   % false for flux units
DefV.MinArea            = 1;
% Secondary parameters for detection
DefV.AreaOpenConn       = 8;  % [4 | 8]
DefV.RegionMaxConn      = 8;  % [4 | 8]
DefV.ReplaceVal         = []; %[-Inf 0];
DefV.ColXY              = {'XPIX_PEAK','YPIX_PEAK'};
% Output columns
DefV.ColFile            = [];
% DefV.ColCell            = {'XWIN_IMAGE','YWIN_IMAGE','XPEAK_IMAGE','YPEAK_IMAGE',...
%                            'X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE',...
%                            'ALPHAWIN_J2000','DELTAWIN_J2000',...
%                            'PEAKF_VAL','PEAKF_VALTOT',...
%                            'SN','SN_UNF','SN_PSF','SN_ADD',...
%                            'DELTA_S',...
%                            'FLUX_PSF','FLUXERR_PSF','MAG_PSF','MAGERR_PSF','PSF_CHI2','PSF_CHI2BACK','PSF_CHI2CR',...
%                            'FLUX_APER','FLUXERR_APER',...
%                          
%                            'BACK','BACK_STD','BACK_ANNULUS','STD_ANNULUS',...
%                            'AZ','ALT','AIRMASS','PARANG',...
%                            'PEAK_GRADBACK',...
%                            'FLAGS',...
%                            'FLUX_APER_IG','FLUXERR_APER_IG',...
%                            'CAREA','CAREA_NSRC',...
%                            'X_IMAGE','Y_IMAGE',...
%                            'X2_IMAGE','Y2_IMAGE','XY_IMAGE',...
%                            'FLUX_ISO','FLUXERR_ISO',...
%                            'NEAREST_SRCIND','NEAREST_SRCDIST'};
                       
DefV.ColCell            = {'ISFORCE','XWIN_IMAGE','YWIN_IMAGE','XPEAK_IMAGE','YPEAK_IMAGE',...
                           'X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE',...
                           'THETA','ELONGATION',...
                           'ALPHAWIN_J2000','DELTAWIN_J2000',...
                           'PEAKF_VAL','PEAKF_VALTOT',...
                           'SN','SN_UNF','SN_PSF','SN_ADD',...
                           'MAG_SN',...
                           'FLUX_PSF','FLUXERR_PSF','MAG_PSF','MAGERR_PSF','PSF_CHI2','PSF_CHI2BACK','PSF_CHI2CR',...
                           'FLUX_APER','FLUXERR_APER',...
                           'BACK','BACK_STD','BACK_ANNULUS','STD_ANNULUS',...
                           'PEAK_GRADBACK',...
                           'FLAGS','NEAREST_SRCDIST'};
                       
% DefV.ColCell            = {'XWIN_IMAGE','YWIN_IMAGE','XPEAK_IMAGE','YPEAK_IMAGE',...
%                            'X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE',...
%                            'THETA','ELONGATION',...
%                            'ALPHAWIN_J2000','DELTAWIN_J2000',...
%                            'PEAKF_VAL','PEAKF_VALTOT',...
%                            'SN','SN_UNF','SN_ADD','FLAGS'};                       
                       
% Force positions
DefV.ForcePos           = [];
DefV.ForceCatCol        = {'XWIN_IMAGE','YWIN_IMAGE'};
DefV.OnlyForce          = false;
DefV.WinPosFromForce    = true;   % windowed position are from ForcePos
% Filters
DefV.FilterFind         = true;
DefV.FilterFindThresh   = 10;
DefV.FilterIm           = true;
DefV.Filter             = []; %@Kernel2.gauss;
DefV.FilterFunPar       = {1.5,1.5,0,15,15};
DefV.DefFilter          = @Kernel2.gauss;
DefV.DefFilterFunPar    = {1.5,1.5,0,15,15};
DefV.FilterFun          = 'auto';
DefV.GetPsfPar          = {};    % additional parameters to pass to getmpsf.m
DefV.AddFilter          = {@Kernel2.exp};    % cell array of filters to use on all images
DefV.AddFilterPar       = {{3,15,15}};
DefV.ConvAddFilter      = true;
% Magnitudes and zero points
DefV.ZP                 = 22;
DefV.MagInLuptitude     = true;
DefV.LuptSoft           = 1e-10;
% PSF photometry
DefV.PsfFitRad          = 2;   % pix
% aperture photometry
DefV.AperRad            = [2 4 8 12 16];
DefV.Annulus            = [20 24];
DefV.AperBackFun        = @median; %@rmean; %@median;
DefV.AperBackFunPar     = {}; %{1,[0.25 0.25]};  %{};
DefV.AperBackErrFun     = @std;
DefV.AperBackErrFunPar  = {};
% Gain
DefV.Gain               = {'GAIN'};  % or numeric
DefV.OrigGainKey        = 'ORIGGAIN';
% Background
DefV.BackSub            = false;   % indicate if back already subtracted
DefV.ReBackStd          = false;
DefV.BackPar            = {'Block',[256 256]};

% PSF estimation
% DefV.PSF_EstimatorPar   = {'NoiseFactor',3,'Rad0',8}; %120,'StampSize',120};
% DefV.PSF_SelectPar      = {'ColSN','SN',...
%                            'ColMaxFlux','PEAKF_VALTOT',...
%                            'ColBitFlag','FLAGS',...
%                            'MaskVal',0,...
%                            'PosCol',{'XWIN_IMAGE','YWIN_IMAGE'},...
%                            'BoundryDist',10,...
%                            'MinSN',15,...
%                            'SatLevel',40000};
DefV.HalfSizePSF          = 10;
DefV.PSFSelectorFun       = @psf_cat_selector;
DefV.PSFSelectorPar       = {'ColSN','SN',...
                           'ColMaxFlux','PEAKF_VALTOT',...
                           'ColBitFlag','FLAGS',...
                           'MaskVal',0,...
                           'PosCol',{'XWIN_IMAGE','YWIN_IMAGE'},...
                           'BoundryDist',10,...
                           'MinSN',10,...
                           'SatLevel',60000.*1.6};
DefV.PSFCombinerPar       = {};

% Properties measurment
DefV.PropThresh         = 1.5;  % [sigma] Threshold from bw connectivity

% Time and position
DefV.JD                 = [];   % if given then will be used
DefV.JuldayPar          = {};
DefV.GeodLon            = 'OBSLON';
DefV.GeodLat            = 'OBSLAT';
DefV.GeodUnits          = 'deg';
DefV.TypeLST            = 'm';     % LST type: {'m'|'a'}



% Source cleaning
DefV.CleanEdge          = true;    % Not applied to forced position sources
DefV.CleanEdgeDist      = 2;       

%DefV.CleanLocalSN       = true;
%DefV.CleanByChi2Back    = true;    % remove sources with: PSF_CHI2BACK<PSF_CHI2 (background)
%DefV.CleanByChi2CR      = true;    % remove sources with: PSF_CHI2CR<PSF_CHI2 unless saturated (CR)
DefV.CleanByBackGrad    = true;    % SizePSF = 3; select by: (SN-Thresh)>PEAK_GRADBACK.*SizePSF./(BACK_STD./SizePSF) - important
DefV.GradBackSize       = 64;

% Moments
DefV.MomAperRad         = 6;   % be careful - decreasing generate artifacts
DefV.MomSigma           = 1.5;
DefV.MomMaxIter         = 5;

% Coordinate units
DefV.OutCooUnits        = 'deg';   % {'deg'|'rad'}

% Nearest neighboor
DefV.NearSearchRad      = 30;   % pix

% Mask image
DefV.Mask               = true;   % true|false|MASK|matrix
DefV.MaskDic            = @MASK.def_bitmask_pipeline;
DefV.FlagClass          = 'uint32';
DefV.FlagRadius         = 3;      % FLAGS are searched within this radius
DefV.FlagOperation      = @Util.array.bitor_array;  % @bitand_array | @bitor_array
DefV.FlagDiffSpike      = true;
DefV.DiffSpikePar       = {};
DefV.CleanDiffSpike     = true;
DefV.Bit_Spike          = 'Bit_Spike';

DefV.SearchCR           = true;
DefV.MethodCR           = {'chi2backcr'}; % {'SN_UNFfit','chi2backcr'};
DefV.CleanCR            = true;
DefV.CleanBitNames      = {'Bit_CR','Bit_CR_MEXTRACTOR'};

% Output
%DefV.OutType            = 'sim';   % 'sim'|'astcat'
DefV.SortOut            = 'YWIN_IMAGE';

% SIM related
DefV.ExecField          = SIM.ImageField;

% Verbose
DefV.Verbose            = true;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% inilalize MetaData (in case not used)
MetaData = [];

SizeSim = size(Sim);
Nsim    = numel(Sim);

%-----------------------------
%--- Read ColFile if given ---
%-----------------------------
if (~isempty(InPar.ColFile))
    % Override 'ColCell'
    [~,InPar.ColCell] = create_list(sprintf('@%s',InPar.ColFile),NaN);
end


%----------------------------
%--- set images gain to 1 ---
%----------------------------
%  Im, BackIm, ErrIm fields
if (InPar.Verbose)
    fprintf('  Set gain to 1\n');
end
Sim = gain_correct(Sim,'GainKeys',InPar.Gain,'ExecField',{ImageField,BackField,ErrField},'ExecFieldSqrt',{},...
                       'CatColUpdate',{},'OrigGainKey',InPar.OrigGainKey);
        


%--------------------------------------------
%--- Check if filter supplied by the user ---
%--------------------------------------------
% If no filter found set to some default filter
if (isempty(InPar.Filter))
    Psf = getmpsf(Sim,InPar.GetPsfPar{:});
    if (any(Util.cell.isempty_cell(Psf)))
        % No Filter and no PSF in SIM
        if (InPar.FilterFind)
            %---------------------------
            %--- Find optimal filter ---
            %---------------------------
            % Execute pre-mextractor run
            % This will be used only in order to look for bright stars
            % for PSF estimation
            if (InPar.Verbose)
               fprintf('mextractor first pass - Constructing optimal PSF\n');
            end
            SimPass1 = mextractor(Sim,'FilterFind',false,...
                           'Thresh',InPar.FilterFindThresh,...
                           'ColCell',{'XWIN_IMAGE','YWIN_IMAGE','SN','SN_UNF',...
                                      'FLAGS','PEAKF_VALTOT','X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE',...
                                      'THETA','ELONGATION',...
                                      'NEAREST_SRCDIST'},...
                           'AddFilter',{},...
                           'CleanByBackGrad',false,...
                           'FlagDiffSpike',false,'CleanDiffSpike',false,...
                           'SearchCR',false,'CleanCR',false,...
                           'BackPar',InPar.BackPar);
            %-----------------------------
            %--- Execute psf_estimator ---
            %-----------------------------
            % Note that 'SelectFunPar' should be a cell array and hence
            % the use of (:)
            
%             [Sim,MetaData.ResPSF] = psf_estimator(SimPass1,'OutType','OrigSim',...
%                                                   'SelectFun',@psf_cat_selector,...
%                                                   InPar.PSF_EstimatorPar{:},...
%                                                   'SelectFunPar',InPar.PSF_SelectPar(:));
%                           
 


            [Sim,ResPSF]=psf_extractor(SimPass1,'StampHalfSize',InPar.HalfSizePSF,...
                                           'PSFSelectorFun',InPar.PSFSelectorFun,...
                                           'PSFSelectorPar',InPar.PSFSelectorPar,...
                                           'PSFCombinerPar',InPar.PSFCombinerPar);
            MetaData.ResPSF = ResPSF;
            
            % Convert to a cell array of PSFs
            %InPar.Filter = getmpsf(Sim);
            Psf = getmpsf(Sim);
            if (InPar.Verbose)
                fprintf('mextractor second pass - Source detection and measurments\n');
             end
        else
            % set to default filter
            if (InPar.Verbose)
                fprintf('  No filter found - Use DefFilter\n');
            end
            %InPar.Filter = InPar.DefFilter(InPar.DefFilterFunPar{:});
            Psf = {InPar.DefFilter(InPar.DefFilterFunPar{:})};
        end
    end
else
    
    Psf = {InPar.Filter(InPar.FilterFunPar{:})};
end
Npsf = numel(Psf);
if ~(Npsf==1 || Npsf==Nsim)
    error('Number of PSF should be 1 or like number of images');
end
        
%-------------------------------------
%--- Set up the additional filters ---
%-------------------------------------
% Additional filters are other filters that will be used for source
% detection, except the default filter (usually the optimal filter for
% point source detection). For example, this may include filters for galaxy
% detection.
Nfilt = numel(InPar.AddFilter) + 1;    
if (~isempty(InPar.AddFilter))
    % The user requested for additional filters
    % get the additional filters
    InPar.AddFilter = Kernel2.cellmat(InPar.AddFilter,InPar.AddFilterPar);
end

    
  
if (InPar.CleanDiffSpike)
    InPar.FlagDiffSpike = true;
end


%-----------------------------------
%--- Read and set the Mask image ---
%-----------------------------------
if (~islogical(InPar.Mask))
    % Convert the mask into a MASK object (array2mask)
    % Copy the MASK object into the Sim object
    Sim = mask2sim(array2mask(Mask,InPar.MaskDic),Sim,InPar.MaskDic);
    % Sim now contains the Mask image
    % set InPar.Mask to true
    InPar.Mask = true;
else
    if (InPar.Mask)
        Sim = mask_init(Sim,InPar.FlagClass);
    end
end





switch lower(InPar.OutCooUnits)
    case 'deg'
        Conv2OutCoo = RAD;
    case 'rad'
        Conv2OutCoo = 1;
    otherwise
        error('Unknown OutCooUnits');
end
        


%---------------------
%--- ForcePos list ---
%---------------------
if (isempty(InPar.ForcePos))
    ForcePos = AstCat;  % empty catalog
else
    if (AstCat.isastcat(InPar.ForcePos))
        % ForcePos is already in AstCat format
        ForcePos = InPar.ForcePos;
    elseif (isnumeric(InPar.ForcePos))
        ForcePos = AstCat;
        ForcePos.(CatField)     = InPar.ForcePos;
        ForcePos.(ColCellField) = InPar.ForceCatCol;
        ForcePos                = colcell2col(ForcePos);
    else
        error('Uknown ForcePos option');
    end
end
Nforce = numel(ForcePos);


%------------------------------------
%--- Calculate background and std ---
%------------------------------------
if (InPar.ReBackStd)
    % User requested to re-calculate background and std
    ExistBack = false(SizeSim);
    ExistStd  = false(SizeSim);
else
    % Use BackIm and ErrIm fields if exist
    ExistBack = isfield_populated(Sim,BackField);
    ExistStd  = isfield_populated(Sim,ErrField);
end

% back/std calculation
if (any(~ExistBack) || any(~ExistStd))
    % some BackIm, ErrIm fields need to be recalculated
    % don't subtract background from image
    if (InPar.Verbose)
        fprintf('  Estimate images background and std\n');
    end
    Sim = background(Sim,InPar.BackPar{:},'SubBack',false);
end

% If background is a scalar set CleanByBackGrad to false
if (InPar.CleanByBackGrad)
    BackSize = imagesize(Sim,true,SIM.BackField);
    if any(prod(BackSize,2)==1)
        InPar.CleanByBackGrad = false;
        warning('CleanByBackGrad was set to false since at least one image has a scalar background');
    end
end
    



%---------------------------------------
%--- Subtract background from images ---
%---------------------------------------
if (~InPar.BackSub)
    % Background was not pre-subtracted by user
    % subtract background
    if (InPar.Verbose)
        fprintf('  Subtract background\n');
    end
    SimB = sub_background(Sim);
else
    % User supplied a background subtracted images
    SimB = Sim;
end
    


%-------------------------------           
%--- Define pre-calculations ---
%-------------------------------
% Define what needed to be calculated given the requested output columns
Calc.CooPeak= any(strcmp(InPar.ColCell,'ALPHAPEAK_J2000')) || any(strcmp(InPar.ColCell,'DELTAPEAK_J2000'));
Calc.CooWin = any(strcmp(InPar.ColCell,'ALPHAWIN_J2000')) || any(strcmp(InPar.ColCell,'DELTAWIN_J2000'));
Calc.AM     = any(strcmp(InPar.ColCell,'AIRMASS'));
Calc.ParAng = any(strcmp(InPar.ColCell,'PARANG'));
Calc.AzAlt  = any(strcmp(InPar.ColCell,'AZ')) || ...
              any(strcmp(InPar.ColCell,'ALT')) || ...
              Calc.AM || Calc.ParAng;
Calc.Mom1   = any(strcmp(InPar.ColCell,'XWIN_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'YWIN_IMAGE')) || ...
              Calc.AzAlt || Calc.CooWin;
Calc.Mom2   = any(strcmp(InPar.ColCell,'X2WIN_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'Y2WIN_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'XYWIN_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'A')) || ...
              any(strcmp(InPar.ColCell,'B')) || ...
              any(strcmp(InPar.ColCell,'THETA')) || ...
              any(strcmp(InPar.ColCell,'ELONGATION')) || ...
              any(strcmp(InPar.ColCell,'ELLIPTICITY'));
Calc.PSF    = any(strcmp(InPar.ColCell,'FLUX_PSF')) || ...
              any(strcmp(InPar.ColCell,'FLUXERR_PSF')) || ...
              any(strcmp(InPar.ColCell,'MAG_PSF')) || ...
              any(strcmp(InPar.ColCell,'MAGERR_PSF')) || ...
              any(strcmp(InPar.ColCell,'PSF_CHI2')) || ...
              any(strcmp(InPar.ColCell,'PSF_CHI2BACK')) || ...
              any(strcmp(InPar.ColCell,'PSF_CHI2CR')) || ...
              any(strcmp(InPar.ColCell,'PSF_BACKERR'));
Calc.Mom1   = Calc.Mom1 || Calc.PSF;    
Calc.AperIL = any(strcmp(InPar.ColCell,'FLUX_APER')) || ...
              any(strcmp(InPar.ColCell,'FLUXERR_APER')) || ...
              any(strcmp(InPar.ColCell,'MAG_APER')) || ...
              any(strcmp(InPar.ColCell,'MAGERR_APER')) || ...
              any(strcmp(InPar.ColCell,'FLUX_APER_IL')) || ...
              any(strcmp(InPar.ColCell,'FLUXERR_APER_IL')) || ...
              any(strcmp(InPar.ColCell,'MAG_APER_IL')) || ...
              any(strcmp(InPar.ColCell,'MAGERR_APER_IL'));
Calc.AperIG = any(strcmp(InPar.ColCell,'FLUX_APER_IG')) || ...
              any(strcmp(InPar.ColCell,'FLUXERR_APER_IG')) || ...
              any(strcmp(InPar.ColCell,'MAG_APER_IG')) || ...
              any(strcmp(InPar.ColCell,'MAGERR_APER_IG'));
% Back and std based on aperture photometry
Calc.AperIb = any(strcmp(InPar.ColCell,'BACK_ANNULUS')) || any(strcmp(InPar.ColCell,'STD_ANNULUS'));

Calc.BWconn = any(strcmp(InPar.ColCell,'CAREA')) || ...
              any(strcmp(InPar.ColCell,'CAREA_NSRC')) || ...
              any(strcmp(InPar.ColCell,'X_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'Y_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'X2_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'Y2_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'XY_IMAGE')) || ...
              any(strcmp(InPar.ColCell,'FLUX_ISO'));
Calc.Near   = any(strcmp(InPar.ColCell,'NEAREST_SRCIND')) || ...
              any(strcmp(InPar.ColCell,'NEAREST_SRCDIST'));
Calc.Flags  = any(strcmp(InPar.ColCell,'FLAGS'));
Calc.BackGrad = any(strcmp(InPar.ColCell,'PEAK_XGRADBACK')) || ...
              any(strcmp(InPar.ColCell,'PEAK_YGRADBACK')) || ...
              any(strcmp(InPar.ColCell,'PEAK_GRADBACK'));
Calc.CleanByBackGrad = any(strcmp(InPar.ColCell,'SN')) && ...
                       any(strcmp(InPar.ColCell,'PEAK_GRADBACK')) && ...
                       any(strcmp(InPar.ColCell,'BACK_STD'));
                       

% Check if CleanByBackGrad can be done:
if (InPar.CleanByBackGrad && ~Calc.CleanByBackGrad)
    warning('CleanByBackGrad was set to false as the user did not requested for the SN, PEAK_GRADBACK, BACK_STD columns');
    InPar.CleanByBackGrad = false;
end

%-------------------------------
%--- Calculate JD for images ---
%-------------------------------
% get time collectivly
if (Calc.AzAlt || Calc.AM || Calc.ParAng)
    % get JD for mid exposure
    if (isempty(InPar.JD))
        % JD is not provided by user
        % get from header
        JD = julday(Sim,InPar.JuldayPar{:});
    else
        % JD is provided by user
        JD = InPar.JD;
    end
    % Convert JD into a row vector
    JD = JD(:).';
    
    % get geodetic position
    GeodPos = geodpos(Sim,'Long',InPar.GeodLon,'Lat',InPar.GeodLat,...
                          'InAngUnits',InPar.GeodUnits);
    % Output:
    % GeodPos.Long - observatory longitude [rad]
    % GeodPos.Lat  - observatory latitude [rad]
    
    % Calculate local sidereal time
    %LST = lst(JD,[GeodPos.Long],InPar.TypeLST);   % [fracday]
    LST = celestial.time.lst(JD,[GeodPos.Long],InPar.TypeLST);   % [fracday] Na'ama, 20180808

end
    

%------------------------------------------------
%--- Expand multi aperture photometry columns ---
%------------------------------------------------
Naper = numel(InPar.AperRad);
if ((Calc.AperIL || Calc.AperIG) && numel(InPar.AperRad)>1)
    
    % duplication of aperture photometry columns for all the requested
    % apertures.
    % Note that here we assume that InPar.ColCell is a row vector
    APER_Cols = {'FLUX_APER','FLUXERR_APER','MAG_APER','MAGERR_APER',...
                 'FLUX_APER_IL','FLUXERR_APER_IL','MAG_APER_IL','MAGERR_APER_IL',...
                 'FLUX_APER_IG','FLUXERR_APER_IG','MAG_APER_IG','MAGERR_APER_IG'};
             
             
    % The following call duplicate some of the columns in ColCell
    % but without adding numbers (e.g., "_1", "_2") to each ColCell name.
    % The addition of numbers is done when the columns are populated
    InPar.ColCell = expand_multiple_cell(InPar.ColCell,APER_Cols,Naper);
    
    % Duplication of SN_ADD (additional filters S/N) columns for all the
    % requested additional filters.
    Filters_Cols = {'SN_ADD'};
    InPar.ColCell = expand_multiple_cell(InPar.ColCell,Filters_Cols,Nfilt-1);
    
    
%     Nac       = numel(APER_Cols);
%     for Iac=1:1:Nac,
%         Ind = find(strcmp(InPar.ColCell,APER_Cols{Iac}));
%         % duplicate the aperture phot columns
%         if (~isempty(Ind)),
%             InPar.ColCell = Util.cell.cell_insert(InPar.ColCell,InPar.ColCell(Ind),Ind.*ones(1,Naper-1));
%         end
%     end
end

Ncol    = numel(InPar.ColCell);    % number of requested columns in the output catalog
      




% allocate the filtered SIM object
SimF  = SIM(Nfilt,1);


%----------------------
%--- for each image ---
%----------------------
Cat = Sim; %AstCat(Nsim,1); % The outpu is always SIM
for Isim=1:1:Nsim
    if (InPar.Verbose)
        fprintf('mextractor image %d out of %d\n',Isim,Nsim);
    end
    
    % Initilaize catalog
    Cat(Isim).(CatField) = [];
    
    % For each SIM element
    
    
    % index for the ForcePos AstCat object
    Iforce = min(Isim,Nforce);
    
    % PSF for current image
    ImagePsf = Psf{min(Isim,Npsf)};

    %----------------------------------------------------------------
    %--- Convolve the additional filters withe the current filter ---
    %----------------------------------------------------------------
    if (~isempty(InPar.AddFilter))    
        if (InPar.ConvAddFilter)
            % The user requested to convolve the additional filters with the
            % default filter (e.g., optimal point source filter).
            % Note that here we assume that InPar.Filter is a 2D array.
            InPar.AddFilter = ImUtil.Im.conv2_cell(InPar.AddFilter,ImagePsf,'auto');
            
        end
    end

    %------------------------
    %--- Filter the image ---
    %------------------------
    % Apply the basic (first) filter:
    if (InPar.Verbose)
        fprintf('  Filter image number %d - use %d filters\n',Isim,Nfilt);
    end
    
    % Creation of the filtered image
    % Note that the filtered images corresponding to SimB(Isim) are stored
    % in SimF - the size of SimF doesn't corresponds to the size of SimB!!
    Ifilt = 1;
    if (InPar.FilterIm)
        % Filter the SIM image
        [SimF(Ifilt)] = filter(SimB(Isim),ImagePsf,'ExecField',InPar.ExecField,...
                                'GetPsfPar',InPar.GetPsfPar,...
                                'FilterFunPar',InPar.FilterFunPar,...
                                'FilterFun',InPar.FilterFun);
    else
        % do not filter the image
        SimF(Ifilt) = SimB(Isim);
        % The PSFmaybe required for PSF photometry
        % Psf = getmpsf(Sim(Isim),InPar.GetPsfPar{:}); --- already in
        % ImagePSF
    end
    % Apply the additional filters
    for Ifilt=2:2:Nfilt
        % filter the image and update the ErrIm field:
        SimF(Ifilt) = filter(SimB(Isim),InPar.AddFilter{Ifilt-1},'ExecField',InPar.ExecField,...
                                'FilterFunPar',InPar.AddFilterPar{Ifilt-1},...
                                'FilterFun',InPar.FilterFun);
    end
    % SimF contains the Sim(Isim) image filtered by the various filters
    
    %--------------------------------------------
    %--- Threshold the image and local maxima ---
    %--------------------------------------------
    % SimF contains the filtered image and the updated ErrIm
    
    [SimFC,SimDet] = threshold(SimF,InPar.Thresh, 'ExecField',    InPar.ExecField,...
                                            'UseErrIm',     InPar.ThreshIsSigma,...
                                            'CombineFilter',true,...
                                            'MinArea',      InPar.MinArea,...
                                            'AreaOpenConn', InPar.AreaOpenConn,...
                                            'RegionMaxConn',InPar.RegionMaxConn,...
                                            'ReplaceVal',   InPar.ReplaceVal,...
                                            'ColXY',        InPar.ColXY);
    
    if (InPar.OnlyForce)
        % use only the user supplied list
        PeakX     = [];
        PeakY     = [];
        PeakForce = [];
    else
        % detect sources also in image
        if (InPar.Verbose)
            fprintf('    Threshold image and locate local maxima\n');
        end
        % call threshold was here - moved before the if statement
        
        ColInd    = colname2ind(SimFC(1),InPar.ColXY);
        
        if (InPar.CleanEdge)
            % Clean sources found near the image edges
            % This is applied only for source found by filtering and not to
            % forced position sources
            SizeIm    = size(Sim(1).(InPar.ExecField));
            Range     = [InPar.CleanEdgeDist, SizeIm(2)-InPar.CleanEdgeDist, InPar.CleanEdgeDist, SizeIm(1)-InPar.CleanEdgeDist];
            SimFC(1)   = select_ccdsec(SimFC(1),Range,ColInd);
        end
        PeakX     = SimFC(1).(CatField)(:,ColInd(1));
        PeakY     = SimFC(1).(CatField)(:,ColInd(2));
        PeakForce = false(size(PeakX));
    end

    % The catalog in SimFC(1) contains the X,Y positions of candidate
    % sources
    
    if (isempty(ForcePos(Iforce).(CatField)))
        % PeakX, PeakY, PeakForce are already populated - do nothing
    else
        ColInd    = colname2ind(ForcePos,InPar.ForceCatCol);
        PeakX     = [PeakX;     ForcePos(Iforce).(CatField)(:,ColInd(1))];
        PeakY     = [PeakY;     ForcePos(Iforce).(CatField)(:,ColInd(2))];
        PeakForce = [PeakForce; true(size(PeakX))];
    end
    
    %--------------------------------------
    %--- Sort all peaks by Y coordinate ---
    %--------------------------------------
    [PeakY,SortInd] = sort(PeakY);
    PeakX           = PeakX(SortInd);
    PeakForce       = PeakForce(SortInd);
    
    
    
    
    % number of sources
    Nsrc = numel(PeakX);
    if (InPar.Verbose)
        fprintf('    Extracting %d sources\n',Nsrc);
    end
    % In this stage PeakX, PeakY, PeakForce are populated
    % Note that PeakX/PeakY based on ForcePos request are no necesserly an
    % integer pixel.
    
    %--------------------------------------------------------
    %--- Perform the measurments on all source candidates ---
    %--------------------------------------------------------
    ColCell = InPar.ColCell;
    
    %---------------------------------------
    %--- function to be called only once ---
    %---------------------------------------
    % note that PeakX/PeakY are rounded to an integer pixel
    %PeaksInd = sub2ind(size(SimB(Isim).(InPar.ExecField)),round(PeakY),round(PeakX));
    SizeSimB = size(SimB(Isim).(InPar.ExecField));
    % A faster implementation of sub2ind:
    PeaksInd = uint32(round(PeakY) + (round(PeakX)-1).*SizeSimB(1));
    
    if (any(PeaksInd<1) || any(PeaksInd>prod(SizeSimB)))
        error('Objects outside image boundries');
    end
    
    
    if (numel(SimB(Isim).(BackField))==1)
        PeaksIndBack = uint32(1);
    else
        PeaksIndBack = PeaksInd;
    end
    if (numel(SimB(Isim).(ErrField))==1)
        PeaksIndErr = uint32(1);
    else
        PeaksIndErr = PeaksInd;
    end
        
    % claculate 1st and 2nd moments
    % for all sources - forced and unforced
    Mom1 = [];
    Mom2 = [];
    if (Calc.Mom2)
        if (InPar.Verbose)
            fprintf('  Calculate 1st and 2nd moments\n');
        end
        [Mom1,Mom2] = ImUtil.Im.im_moments(SimB(Isim).(InPar.ExecField),PeakX,PeakY,InPar.MomAperRad,InPar.MomSigma,InPar.MomMaxIter);
    else
        if (Calc.Mom1 || Calc.CooWin || Calc.Near)
            if (InPar.Verbose)
                fprintf('    Calculate 1st moments\n');
            end
            [Mom1]      = ImUtil.Im.im_moments(SimB(Isim).(InPar.ExecField),PeakX,PeakY,InPar.MomAperRad,InPar.MomSigma);
        else
            % Moments were not requested
        end
    end
    % use user forced position instead of Mom1.X, Mom1.Y
    if (InPar.WinPosFromForce)
        Mom1.X(PeakForce) = PeakX(PeakForce);
        Mom1.Y(PeakForce) = PeakY(PeakForce);
    end
    
    % PSF photometry
    if (Calc.PSF)
        % The position of forced objects are taken from the forced position
        if (InPar.Verbose)
            fprintf('    Calculate PSF photometry\n');
        end
        
        % Fit only the PSF core
        % Cut the PSF core:
        
        % should change psf_phot to SIM/psf_phot
%         [PsfPhot,PsfPhotCol] = psf_phot(Sim(Isim).(InPar.ExecField),Sim(Isim).(BackField),ImagePsf,[Mom1.X,Mom1.Y],...
%                                     'FitRad',InPar.PsfFitRad,...
%                                     'NormPSF',true,...
%                                     'OutType','mat');
        [PsfPhot,PsfPhotCol] = psf_phot(Sim(Isim),[Mom1.X,Mom1.Y],...
                                    'FitRad',InPar.PsfFitRad,...
                                    'NormPSF',true);
                                
                 
    end
    
    % aperture photometry (IL)
    if (Calc.AperIL || Calc.AperIb)
        % The position of forced objects are taken from the forced
        % positions
        if (InPar.Verbose)
            fprintf('    Calculate IL aperture photometry (using aper_phot)\n');
        end
%         [AperPhotIL,AperPhotCol] = aper_phot(Sim(Isim).(InPar.ExecField),[Mom1.X,Mom1.Y],...
%                                                 'AperRad',InPar.AperRad,...
%                                                 'Annulus',InPar.Annulus,...
%                                                 'BackFun',InPar.AperBackFun,...
%                                                 'BackFunPar',InPar.AperBackFunPar,...
%                                                 'BackErrFun',InPar.AperBackErrFun,...
%                                                 'BackErrFunPar',InPar.AperBackErrFunPar,...
%                                                 'Back',[]);
        [AperPhotIL,AperPhotCol] = aper_phot(Sim(Isim),[Mom1.X,Mom1.Y],...
                                                'AperRad',InPar.AperRad,...
                                                'Annulus',InPar.Annulus,...
                                                'BackFun',InPar.AperBackFun,...
                                                'BackFunPar',InPar.AperBackFunPar,...
                                                'BackErrFun',InPar.AperBackErrFun,...
                                                'BackErrFunPar',InPar.AperBackErrFunPar,...
                                                'Back',[]);
        
    end
    
    % aperture photometry (IG)
    if (Calc.AperIG)
        % The position of forced objects are taken from the forced
        % positions
        if (InPar.Verbose)
            fprintf('    Calculate IG aperture photometry (using aper_phot)\n');
        end
        [AperPhotIG,AperPhotCol] = aper_phot(Sim(Isim),[Mom1.X,Mom1.Y],...
                                                'AperRad',InPar.AperRad,...
                                                'Annulus',InPar.Annulus,...
                                                'BackFun',InPar.AperBackFun,...
                                                'BackFunPar',InPar.AperBackFunPar,...
                                                'BackErrFun',InPar.AperBackErrFun,...
                                                'BackErrFunPar',InPar.AperBackErrFunPar,...
                                                'Back',SimB(Isim).(BackField)(PeaksIndBack));
        
        
    end
    
    % WCS positions
    % % fixing a bug in WCS - need to fix this within SIM...
    Wsim = ClassWCS.populate(Sim(Isim));
   
    if (Calc.CooPeak)
        [AlphaPeak,DeltaPeak] = xy2coo(Wsim,[PeakX,PeakY]);   % radians
        
        %[AlphaPeak,DeltaPeak] = xy2coo(Sim(Isim),PeakX,PeakY);   % radians
    end
    if (Calc.CooWin)
        [AlphaWin,DeltaWin] = xy2coo(Wsim,[Mom1.X,Mom1.Y]);   % radians
        
        %[AlphaWin,DeltaWin] = xy2coo(Sim(Isim),Mom1.X,Mom1.Y);   % radians
    end
        
    
    % get time - this was done collectivly for all the images outside the
    % loop.
    if (Calc.AzAlt || Calc.AM || Calc.ParAng)
        % calculate local sidereal time
        
        % added by Na'ama, 20180808
        if (Calc.CooWin)
            Alpha = AlphaWin;
            Delta = DeltaWin;
        else
            Alpha = AlphaPeak;
            Delta = DeltaPeak;
        end
        
        
        % DeltaWin, AlphaWin [radians]
        % LST [fracday]
        [Az,Alt,AirMass] = celestial.coo.ha2az(LST(Isim).*2.*pi-Alpha,Delta,GeodPos(Isim).Lat);   % [rad,rad,[]] Na'ama, 20180808
        
        if (Calc.ParAng)
            ParAng           = celestial.coo.parallactic_angle(Alpha, Delta, LST(Isim), GeodPos(Isim).Lat);  % [rad] Na'ama, 20180808
        end
    end
   
    if (Calc.BWconn)
        % Threshold the image above the InPar.PropThresh parameter
        DetConn = SimB(Isim).(InPar.ExecField)./SimB(Isim).(ErrField) > InPar.PropThresh;
        % Generate pixel connectivity lists
        % These are list of pixels of connected regions above the
        % InPar.PropThresh parameter
        CC = bwconncomp(DetConn,4);
     
        % Connectivity labeled image
        LabeledImage = labelmatrix(CC);
        
        % These sources are not related to the sources found
        % by the thresholding - so first we need to associate them with the
        % sources found earlier
        
        % Get the label index in CC for each peak index in PeaksInd
        PeaksNum     = (1:1:Nsrc).';           % serial number of peak
        IndexInCC    = LabeledImage(PeaksInd); % index of peak in CC
        Flag0        = IndexInCC==0;           % Peaks with no ID label matrix
        IndexInCC    = IndexInCC(~Flag0);      % index of peaks in CC (only with ID)
        IndexInPeaks = PeaksNum(~Flag0);       % index of peaks in PeaksInd (only with ID)
        
        % Implement the following in the prop calc case  block
        Nsrc_prop   = numel(IndexInCC);
        SizeIm      = imagesize(Sim(Isim));
        [MatX,MatY] = meshgrid((1:1:SizeIm(1)),(1:1:SizeIm(2)));
        
    end
    
    
    if (Calc.Near)
        [XY,SI] = sortrows([Mom1.X, Mom1.Y],2);
        % need a better search_cat.m (nicer for the user)
        ResNear = VO.search.search_cat(XY,XY,[],'CooType','plane',...
                                                        'SearchMethod','binmdup',...
                                                        'SearchRad',InPar.NearSearchRad);
    
    end
    
    %-----------------------------------
    %--- Diffraction spikes flagging ---
    %-----------------------------------
    if (InPar.FlagDiffSpike)
        % Populate the SIM catalog with temporary X,Y coordinates
        % This is used only for the spikes finding
        Sim(Isim).(CatField) = [PeakX,PeakY];
        Sim(Isim).(ColCellField) = {'X','Y'};
        Sim(Isim) = colcell2col(Sim(Isim));
        % spikes finding
        [Sim(Isim),~] = find_diff_spikes(Sim(Isim),InPar.DiffSpikePar{:});
        % Copy the MASK image to the output Cat object
        Cat(Isim).(MaskField) = Sim(Isim).(MaskField);
        
    end

    
    
    %---------------------------------------
    %--- Calculate Flags for each source ---
    %---------------------------------------
    if (Calc.Flags)     
        Flags = zeros(Nsrc,1,class(Sim(Isim).(MaskField)));
    end
    
%     if (Calc.Flags)           
%         if (InPar.Mask && isfield_populated(Sim(Isim),MaskField))
%             % If Mask is true then propogate the Mask into the source Flags
%             % find the indices of pixels that belongs to each source
%             CellSrcInd = ImUtil.Im.find_within_radius_cell(SizeSimB,PeakX,PeakY,InPar.FlagRadius,true);
%             % Calculate the integrated Flag (within FlagRadius) for each
%             % source
%             Flags = mask4src(Sim(Isim),CellSrcInd,InPar.FlagOperation);
%         else
%             % do not propogate the Mask
%             Flags = zeros(Nsrc,1,InPar.FlagClass);
%         end
%     end
    
    
    %---------------------------------------------------------
    %--- Calculate the background gradient for each source ---
    %---------------------------------------------------------
    if (Calc.BackGrad)
        if (numel(SimB(Isim).(BackField))>1)
            % Background image is an image
            [BackGradX,BackGradY] = gradient(SimB(Isim).(BackField));
        else
            % Background image is a scalar
            % can not calculate gradient
            % populate with NaNs.
            BackGradX = NaN;
            BackGradY = NaN;
        end
    end
    
    
    
    %------------------------------------
    %--- Populate the catalog columns ---
    %------------------------------------
    % Cat(Isim).(CatField) = nan(Nsrc,Ncol); - Ncol is unknown yet??
    Icol = 0;
    %for Icol=1:1:Ncol,
    while Icol<Ncol
        Icol = Icol + 1;
        % HOW TO ADD A NEW PROPERTY
        %---------------------------
        % In order to calculate a property you have to:
        % 0. InArg.unlearn('mextractor');
        % 1. define it in InPar.ColCell
        % 2. Add the treatment to the specific ColCell to the switch/case
        % below - (i.e., Cat(Isim).(CatField)(:,Icol) = Property_Values)
        % Note that Cat is an AstCat object containing the output catalog,
        % Isim is the image index, CatField is a string indicating the
        % Catalog field name in the AstCat object and Icol is the index of
        % the column you want to populate in the Cat object.
        % 3. Property_Values must be either a scalar or a column vector of  
        % the same length as the lines in the catalog.
        % 4. You can either calculate Property_Values within the "case" or
        % the "function to be called only once" block.
        switch (ColCell{Icol})
            %----------------------------
            %--- Auxilary information ---
            %----------------------------
            case 'NUMBER'
                % serial index of source
                Cat(Isim).(CatField)(:,Icol) = (1:1:Nsrc)';
            case 'ISFORCE'
                % A logical flag indicating the source origin
                % ForcePos list (true).
                Cat(Isim).(CatField)(:,Icol) = PeakForce;
            
            %--------------------------------------------
            %--- Source position in IMAGE coordinates ---
            %--------------------------------------------
            case 'XPEAK_IMAGE'
                % The whole pixel value position
                Cat(Isim).(CatField)(:,Icol) = PeakX;
            case 'YPEAK_IMAGE'
                % The whole pixel value position
                Cat(Isim).(CatField)(:,Icol) = PeakY;
            case 'XWIN_IMAGE'
                % 1st moment with Gaussian weight function
                Cat(Isim).(CatField)(:,Icol) = Mom1.X;
            case 'YWIN_IMAGE'
                % 1st moment with Gaussian weight function
                Cat(Isim).(CatField)(:,Icol) = Mom1.Y;
    
            %-----------------------------
            %--- Source second moments ---
            %-----------------------------
            case 'X2WIN_IMAGE'
                % 2nd moment with Gaussian weight function
                Cat(Isim).(CatField)(:,Icol) = Mom2.X2;
            case 'Y2WIN_IMAGE'
                % 2nd moment with Gaussian weight function
                Cat(Isim).(CatField)(:,Icol) = Mom2.Y2;    
            case 'XYWIN_IMAGE'
                % 2nd moment with Gaussian weight function
                Cat(Isim).(CatField)(:,Icol) = Mom2.XY;    
    
            %-----------------------------------------
            %--- Shape parameters from 2nd moments ---
            %-----------------------------------------
            case 'A'
                % Semi major axis [pix]
                Cat(Isim).(CatField)(:,Icol) = sqrt(0.5.*(Mom2.X2+Mom2.Y2) + sqrt( (Mom2.X2-Mom2.Y2).^2./4 + Mom2.XY.^2 ));
            case 'B'
                % Semi minor axis [pix]
                B = sqrt(0.5.*(Mom2.X2+Mom2.Y2) - sqrt( (Mom2.X2-Mom2.Y2).^2./4 + Mom2.XY.^2 ));
                B(imag(B)~=0) = NaN;
                Cat(Isim).(CatField)(:,Icol) = B;
                
            case 'THETA'
                % Orientation [deg] (-45..45)
                Cat(Isim).(CatField)(:,Icol) =    0.5.*atan(2.*Mom2.XY./(Mom2.X2 - Mom2.Y2)).*RAD;
                
            case 'ELONGATION'
                % Elongation (A/B)
                A = sqrt(0.5.*(Mom2.X2+Mom2.Y2) + sqrt( (Mom2.X2-Mom2.Y2).^2./4 + Mom2.XY.^2 ));
                B = sqrt(0.5.*(Mom2.X2+Mom2.Y2) - sqrt( (Mom2.X2-Mom2.Y2).^2./4 + Mom2.XY.^2 ));
                B(imag(B)~=0) = NaN;
                Cat(Isim).(CatField)(:,Icol) = real(A./B);
                
            case 'ELLIPTICITY'
                % Ellipticity (1-B/A)
                A = sqrt(0.5.*(Mom2.X2+Mom2.Y2) + sqrt( (Mom2.X2-Mom2.Y2).^2./4 + Mom2.XY.^2 ));
                B = sqrt(0.5.*(Mom2.X2+Mom2.Y2) - sqrt( (Mom2.X2-Mom2.Y2).^2./4 + Mom2.XY.^2 ));
                B(imag(B)~=0) = NaN;
                Cat(Isim).(CatField)(:,Icol) = 1 - real(B./A);
            
                
            %-----------------------
            %--- Signal-to-noise ---
            %-----------------------
            case 'SN'
                % S/N in the primary filtered image
                % Note tht SimFC(1) contains the image filtered with the
                % Filter - the other elements of SimF are for the
                % additional filters corresponding to image Isim.
                Cat(Isim).(CatField)(:,Icol) = SimFC(1).(InPar.ExecField)(PeaksInd)./ ...
                                               SimFC(1).(ErrField)(PeaksIndErr);
            case 'MAG_SN'
                % -2.5.*log10(SN)
                % This is similar to PSF magnitude (without ZP)
                Cat(Isim).(CatField)(:,Icol) = -2.5.*log10(SimFC(1).(InPar.ExecField)(PeaksInd)./ ...
                                                           SimFC(1).(ErrField)(PeaksIndErr));
            case 'SN_UNF'
                % S/N in the unfiltered image
                Cat(Isim).(CatField)(:,Icol) = SimB(Isim).(InPar.ExecField)(PeaksInd)./ ...
                                               SimB(Isim).(ErrField)(PeaksIndErr);
            case 'SN_ADD'
                % S/N in the additional filtered images
                if (Ifilt==1)
                    % the user didn't speficy any additional filters
                    % populated SN_ADD with NaN
                    Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                else
                    % The user specified additional filters
                    % Calculate the detection S/N for each of the
                    % additional filters
                    Icol = Icol - 1;
                    for Ifilt=2:1:Nfilt
                        Icol = Icol + 1;
                        Cat(Isim).(CatField)(:,Icol) = SimFC(Ifilt).(InPar.ExecField)(PeaksInd)./ ...
                                                       SimFC(Ifilt).(ErrField)(PeaksIndErr);
                    end
                    
                    % update the ColCell
                    ColCell(Icol:Icol+Nfilt-2) = add_index2cell(ColCell(Icol:Icol+Nfilt-2));
                    % No need to skip columns as Icol was already increased during
                    % the loop.
                    
                    
                end
            
            case 'DELTA_S'
                % Delta score between PSF hypothesis and first hypothesis
                Cat(Isim).(CatField)(:,Icol) = SimF(1).(InPar.ExecField)(PeaksInd) - SimF(2).(InPar.ExecField)(PeaksInd);
                
                K_G = padarray(InPar.AddFilter{Ifilt-1},[3 3],'both')-ImagePsf;
                S2  = filter2(K_G.^2,Sim(Isim).(ErrField).^2,'same');
                %Cat(Isim).(CatField)(:,Icol) = S2(PeaksInd)
            %-----------------------
            %--- Peak properties ---
            %-----------------------
            case 'PEAKF_VAL'
                % Maximum flux at peak pixel (background subtrcated)
                Cat(Isim).(CatField)(:,Icol) = SimB(Isim).(InPar.ExecField)(PeaksInd);
            case 'PEAKF_VALTOT'
                % Maximum flux at peak pixel (including background)
                Cat(Isim).(CatField)(:,Icol) = Sim(Isim).(InPar.ExecField)(PeaksInd);
            case 'MAG_PEAK'
                % Maximum mag at peak pixel
                Cat(Isim).(CatField)(:,Icol) = convert.flux2mag(SimB(Isim).(InPar.ExecField)(PeaksInd),...
                                                        InPar.ZP, InPar.MagInLuptitude,InPar.LuptSoft);
              
            case 'PEAK_XGRAD'
                % Gradient at peak
            case 'PEAK_YGRAD'
                % Gradient at peak
            case 'PEAK_GRAD'
                % Gradient at peak
            case 'PEAK_DEL2'
                % Laplacian at peak
                
            case 'PEAK_XGRADFILT'
                % Gradient at peak of filtered image
            case 'PEAK_YGRADFILT'
                % Gradient at peak of filtered image 
            case 'PEAK_GRADFILT'
                % Gradient at peak of filtered image 
            case 'PEAK_DEL2FILT'
                % Laplacian at peak of filtered image
                
            case 'PEAK_XGRADBACK'
                % Gradient at peak of background image
                Cat(Isim).(CatField)(:,Icol) = BackGradX(PeaksIndBack);
            case 'PEAK_YGRADBACK'
                % Gradient at peak of background image 
                Cat(Isim).(CatField)(:,Icol) = BackGradY(PeaksIndBack);
            case 'PEAK_GRADBACK'
                % Gradient at peak of background image
                Cat(Isim).(CatField)(:,Icol) = sqrt(BackGradX(PeaksIndBack).^2 + BackGradY(PeaksIndBack).^2);
            case 'PEAK_DEL2BACK'
                % Laplacian at peak of background image
                
            %---------------------------
            %--- Aperture photometry ---
            %---------------------------
            % There are several types of aperture photometry:
            % FLUX_APER || FLUX_APER_IL - aperture phortometry based on
            %                             whole pixels (aper_phot.m) with
            %                             background estimated from local
            %                             annulus.
            % FLUX_APER_IG              - aperture phortometry based on
            %                             whole pixels (aper_phot.m) with
            %                             background estimated from global
            %                             background.
            % FLUX_APER_CL              - Convolution based aperture
            %                             photometry with background
            %                             estimated from local annulus.
            % FLUX_APER_CG              - Convolution based aperture
            %                             photometry with background
            %                             estimated from global background.
            case {'FLUX_APER','FLUX_APER_IL'}
                % aperture photometry
                % treat multiple apertures
                Cat(Isim).(CatField)(:,Icol:Icol+Naper-1) = AperPhotIL.(CatField)(:,AperPhotCol.Aper);
                ColCell(Icol:Icol+Naper-1) = add_index2cell(ColCell(Icol:Icol+Naper-1));
                % skip columns - already populated
                Icol = Icol + Naper - 1;
            case {'FLUXERR_APER','FLUXERR_APER_IL'}
                % aperture photometry
                % error is taken as sqrt of flux (assume gain=1)
                Cat(Isim).(CatField)(:,Icol:Icol+Naper-1) = AperPhotIL.(CatField)(:,AperPhotCol.AperErr);
                ColCell(Icol:Icol+Naper-1) = add_index2cell(ColCell(Icol:Icol+Naper-1));
                % skip columns - already populated
                Icol = Icol + Naper - 1;
            case {'MAG_APER','MAG_APER_IL'}
                % aperture photometry
                Cat(Isim).(CatField)(:,Icol:Icol+Naper-1) = convert.flux2mag(AperPhotIL.(CatField)(:,AperPhotCol.Aper),...
                                                        InPar.ZP, InPar.MagInLuptitude,InPar.LuptSoft);
                ColCell(Icol:Icol+Naper-1) = add_index2cell(ColCell(Icol:Icol+Naper-1));
                % skip columns - already populated
                Icol = Icol + Naper - 1;
                
            case {'MAGERR_APER','MAGERR_APER_IL'}
                % aperture photometry
                Cat(Isim).(CatField)(:,Icol:Icol+Naepr-1) = SN2DM .* AperPhotIL.(CatField)(:,AperPhotCol.AperErr)./...
                                                                     AperPhotIL.(CatField)(:,AperPhotCol.Aper);
                ColCell(Icol:Icol+Naper-1) = add_index2cell(ColCell(Icol:Icol+Naper-1));
                % skip columns - already populated
                Icol = Icol + Naper - 1;
                
            case 'BACK_ANNULUS'
                % background in annulus
                Cat(Isim).(CatField)(:,Icol) = AperPhotIL.(CatField)(:,AperPhotCol.Back);
            case 'STD_ANNULUS'
                % std in annulus
                Cat(Isim).(CatField)(:,Icol) = AperPhotIL.(CatField)(:,AperPhotCol.BackErr);
                
            case {'FLUX_APER_IG'}
                % aperture photometry
                % treat multiple apertures
                Cat(Isim).(CatField)(:,Icol:Icol+Naper-1) = AperPhotIG.(CatField)(:,AperPhotCol.Aper);
                ColCell(Icol:Icol+Naper-1) = add_index2cell(ColCell(Icol:Icol+Naper-1));
                % skip columns - already populated
                Icol = Icol + Naper - 1;
            case {'FLUXERR_APER_IG'}
                % aperture photometry
                % error is taken as sqrt of flux (assume gain=1)
                Cat(Isim).(CatField)(:,Icol:Icol+Naper-1) = AperPhotIG.(CatField)(:,AperPhotCol.AperErr);
                ColCell(Icol:Icol+Naper-1) = add_index2cell(ColCell(Icol:Icol+Naper-1));
                % skip columns - already populated
                Icol = Icol + Naper - 1;
            case {'MAG_APER_IG'}
                % aperture photometry
                Cat(Isim).(CatField)(:,Icol:Icol+Naper-1) = convert.flux2mag(AperPhotIG.(CatField)(:,AperPhotCol.Aper),...
                                                        InPar.ZP, InPar.MagInLuptitude,InPar.LuptSoft);
                ColCell(Icol:Icol+Naper-1) = add_index2cell(ColCell(Icol:Icol+Naper-1));
                % skip columns - already populated
                Icol = Icol + Naper - 1;
                
            case {'MAGERR_APER_IG'}
                % aperture photometry
                Cat(Isim).(CatField)(:,Icol:Icol+Naepr-1) = SN2DM .* AperPhotIG.(CatField)(:,AperPhotCol.AperErr)./AperPhotIG.(CatField)(:,AperPhotCol.Aper);
                ColCell(Icol:Icol+Naper-1) = add_index2cell(ColCell(Icol:Icol+Naper-1));
                % skip columns - already populated
                Icol = Icol + Naper - 1;    
            %----------------------------------                                        
            %--- Background level and noise ---
            %----------------------------------
            % additional background are available within aperture phot
            case 'BACK'
                % The background level from the 'BackIm' field
                Cat(Isim).(CatField)(:,Icol) = SimB(Isim).(BackField)(PeaksIndBack);
            case 'BACK_STD'
                % The background noise from the 'ErrIm' field
                Cat(Isim).(CatField)(:,Icol) = SimB(Isim).(ErrField)(PeaksIndErr);
          
            
                
            %-----------------------------------
            %--- Connected pixels properties ---
            %-----------------------------------
            % use bwconncomp.m and regionprops.m
            case 'CAREA'
                % Connected area [pixels]
                % Number of connected pixels above InPar.PropThresh that
                % contains the source.
                % NaN - if a connected region associated with the source was not
                % found
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                Cat(Isim).(CatField)(IndexInPeaks,Icol) = cellfun(@numel,CC.PixelIdxList(IndexInCC));

                warning('NEED to test CAREA option - plots...')
                
            case 'CAREA_NSRC'
                % sources which their connected area cover another source.
                % Number of sources that are enclosed by this connected
                % area.
                % NaN - if source doesnot belong to a connected region.
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                for Isrc_prop=1:1:Nsrc_prop
                    Cat(Isim).(CatField)(IndexInPeaks(Isrc_prop),Icol) = numel(find(IndexInCC==IndexInCC(Isrc_prop)));
                end
                
                warning('NEED to test CAREA_NSRC option - plots... faster using sum(Flag)?')
        
            case 'X_IMAGE'
                % first moment in X (weighted by background subtrcated flux)
                % calculated in the connected area
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                for Isrc_prop=1:1:Nsrc_prop
                    Cat(Isim).(CatField)(IndexInPeaks(Isrc_prop),Icol) = ...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                                    MatX(CC.PixelIdxList{IndexInCC(Isrc_prop)}))./...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}));
                end

            case 'Y_IMAGE'
                % first moment in Y (weighted by background subtrcated flux)
                % calculated in the connected area
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                for Isrc_prop=1:1:Nsrc_prop
                    Cat(Isim).(CatField)(IndexInPeaks(Isrc_prop),Icol) = ...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                                    MatY(CC.PixelIdxList{IndexInCC(Isrc_prop)}))./...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}));
                end    
           
            case 'X2_IMAGE'
                % second moment in X (weighted by background subtrcated flux)
                % calculated in the connected area
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                for Isrc_prop=1:1:Nsrc_prop
                    X1 = sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                             MatX(CC.PixelIdxList{IndexInCC(Isrc_prop)}))./...
                         sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}));
                            
                    Cat(Isim).(CatField)(IndexInPeaks(Isrc_prop),Icol) = ...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                                    MatX(CC.PixelIdxList{IndexInCC(Isrc_prop)}).^2)./...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)})) - X1.^2;
                end
            case 'Y2_IMAGE'
                % second moment in Y (weighted by background subtrcated flux)
                % calculated in the connected area
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                for Isrc_prop=1:1:Nsrc_prop
                    Y1 = sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                             MatY(CC.PixelIdxList{IndexInCC(Isrc_prop)}))./...
                         sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}));
                          
                    Cat(Isim).(CatField)(IndexInPeaks(Isrc_prop),Icol) = ...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                                    MatY(CC.PixelIdxList{IndexInCC(Isrc_prop)}).^2)./...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)})) - Y1.^2;
                end
            case 'XY_IMAGE'
                % second moment in X*Y (weighted by background subtrcated flux)
                % calculated in the connected area
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                for Isrc_prop=1:1:Nsrc_prop
                    X1 = sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                             MatX(CC.PixelIdxList{IndexInCC(Isrc_prop)}))./...
                         sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}));
                    Y1 = sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                             MatY(CC.PixelIdxList{IndexInCC(Isrc_prop)}))./...
                         sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}));
                                
                    Cat(Isim).(CatField)(IndexInPeaks(Isrc_prop),Icol) = ...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).* ...
                                    MatX(CC.PixelIdxList{IndexInCC(Isrc_prop)}).*MatY(CC.PixelIdxList{IndexInCC(Isrc_prop)})   )./...
                                sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)})) - X1.*Y1;
                end
                
            
            % Iso photometry
            case 'FLUX_ISO'
                % Backggropund subtracted flux within connected area
                % NaN - if a source is not associated with a connected area
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                for Isrc_prop=1:1:Nsrc_prop
                    Cat(Isim).(CatField)(IndexInPeaks(Isrc_prop),Icol) = ...
                            sum(SimB(Isim).(InPar.ExecField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}));
                end
                
            case 'FLUXERR_ISO'
                % flux error within connected area
                % NaN - if a source is not associated with a connected area
                % Note that this is calculated assuming the noise is
                % Poisson.
                Cat(Isim).(CatField)(:,Icol) = nan(Nsrc,1);
                for Isrc_prop=1:1:Nsrc_prop
                    Cat(Isim).(CatField)(IndexInPeaks(Isrc_prop),Icol) = ...
                            sqrt(sum(Sim(Isim).(ErrField)(CC.PixelIdxList{IndexInCC(Isrc_prop)}).^2));
                end
                
                
            
                
                
            %----------------------
            %--- PSF photometry ---
            %----------------------
            % Note that the fields of PsfPhotCal:
            % .Flux, .FluxErr, .BackErr, .Chi2 are hard coded in psf_phot.m
            case 'SN_PSF'
                % detection S/N of PSF photometry
                Cat(Isim).(CatField)(:,Icol) = PsfPhot.(CatField)(:,PsfPhotCol.Flux)./AperPhotIL.(CatField)(:,AperPhotCol.BackErr);
                
            case 'FLUX_PSF'
                % PSF photometry
                Cat(Isim).(CatField)(:,Icol) = PsfPhot.(CatField)(:,PsfPhotCol.Flux);
            case 'FLUXERR_PSF'
                % PSF photometry
                Cat(Isim).(CatField)(:,Icol) = PsfPhot.(CatField)(:,PsfPhotCol.FluxErr);
            case 'MAG_PSF'
                Cat(Isim).(CatField)(:,Icol) = convert.flux2mag(PsfPhot.(CatField)(:,PsfPhotCol.Flux),...
                                                        InPar.ZP, InPar.MagInLuptitude,InPar.LuptSoft);
            case 'MAGERR_PSF'
                
                Cat(Isim).(CatField)(:,Icol) = SN2DM .* PsfPhot.(CatField)(:,PsfPhotCol.FluxErr)./PsfPhot.(CatField)(:,PsfPhotCol.Flux);
            case 'PSF_CHI2'
                % PSF photometry \chi2 for point source
                Cat(Isim).(CatField)(:,Icol) = PsfPhot.(CatField)(:,PsfPhotCol.Chi2);
            case 'PSF_CHI2BACK'
                % PSF photometry \chi2 for no source
                Cat(Isim).(CatField)(:,Icol) = PsfPhot.(CatField)(:,PsfPhotCol.Chi2Back);
            case 'PSF_CHI2CR'
                % PSF photometry \chi2 for CR
                Cat(Isim).(CatField)(:,Icol) = PsfPhot.(CatField)(:,PsfPhotCol.Chi2CR);    
            case 'PSF_BACKERR'
                % PSF photometry
                Cat(Isim).(CatField)(:,Icol) = PsfPhot.(CatField)(:,PsfPhotCol.BackErr);
                
            
            %-----------------------------
            %--- WCS J2000 coordinates ---
            %-----------------------------
            case 'ALPHAPEAK_J2000'
                % WCS position for peak
                Cat(Isim).(CatField)(:,Icol) = AlphaPeak.*Conv2OutCoo;
            case 'DELTAPEAK_J2000'
                % WCS position for peak
                Cat(Isim).(CatField)(:,Icol) = DeltaPeak.*Conv2OutCoo;
            case 'ALPHAWIN_J2000'
                % WCS position for Windowed position
                Cat(Isim).(CatField)(:,Icol) = AlphaWin.*Conv2OutCoo;
            case 'DELTAWIN_J2000'
                % WCS position for Windowed position
                Cat(Isim).(CatField)(:,Icol) = DeltaWin.*Conv2OutCoo;
                
                
            %--------------------------------------------------------
            %--- Horizontal coordinates, telescope and atmosphere ---
            %--------------------------------------------------------
            case 'AZ'
                % Azimuth [same units as ALPHA/DELTA]
                Cat(Isim).(CatField)(:,Icol) = Az.*Conv2OutCoo;
                
            case 'ALT'
                % Altitude [same units as ALPHA/DELTA]
                Cat(Isim).(CatField)(:,Icol) = Alt.*Conv2OutCoo;
                
            case 'AIRMASS'
                % Airmass
                Cat(Isim).(CatField)(:,Icol) = AirMass;
                
            case 'PARANG'
                % parallactic angle [same units as ALPHA/DELTA]
                Cat(Isim).(CatField)(:,Icol) = ParAng.*Conv2OutCoo;
                
            %-------------
            %--- Flags ---
            %-------------
            case 'FLAGS'
                % Flags
                % In this stage Flags contain zero
                % populated at the end
                Cat(Isim).(CatField)(:,Icol) = Flags;
                
            %----------------------
            %--- Nearby sources ---
            %----------------------
            case 'NEAREST_SRCIND'
                % Index of nearest source (i.e., Number of source prior to
                % sort) - The index of source is stored in the 'NUMBER'
                % column.
                % NaN if no source within search radius (InPar.NearSearchRad)
                Cat(Isim).(CatField)(:,Icol) = NaN;
                FlagNfound = [ResNear.Nfound]>0;
                Cat(Isim).(CatField)(FlagNfound,Icol) = cat(1,ResNear.IndCat); 
            case 'NEAREST_SRCDIST'
                % Distance to nearest source [pixels]
                % NaN if no source within search radius (InPar.NearSearchRad)
                Cat(Isim).(CatField)(:,Icol) = NaN;
                FlagNfound = [ResNear.Nfound]>0; 
                Cat(Isim).(CatField)(FlagNfound,Icol) = cat(1,ResNear.DistRAD);   
                
                if (any(cat(1,ResNear.DistRAD))<1)
                    warning('There are some distance smaller than 1 - check why');
                end
                  
            otherwise
                error('Unknown requested column: %s in the ColCell input argument',ColCell{Icol});
        end
    end
    % populate the CellCol information
    Cat(Isim).(ColCellField) = ColCell;
    Cat(Isim)                = colcell2col(Cat(Isim));
    
 
    
    
    %------------------------------------
    %--- Post search catalog cleaning ---
    %------------------------------------
    
%     if (InPar.CleanLocalSN),
%         LocalSN = col_get(Cat(Isim),'FLUX_APER_1')./col_get(Cat(Isim),'STD_ANNULUS');
%         Cat(Isim) = row_select(Cat(Isim),LocalSN>InPar.Thresh);
%     end
        
    if (InPar.CleanByBackGrad)
        %SizePSF      = curve_growth_psf(
        
        %SizePSF      = 3;
        %FlagBackGrad = (col_get(Cat(Isim),'SN')-InPar.Thresh)>(col_get(Cat(Isim),'PEAK_GRADBACK').*SizePSF./(col_get(Cat(Isim),'BACK_STD')./SizePSF));
        %Cat(Isim)    = row_select(Cat(Isim),FlagBackGrad);
          
        % Important: current version works only when units are electrons
        BB = col_get(Cat(Isim),'BACK');
        BS = col_get(Cat(Isim),'BACK_STD');
        FlagBackGrad = (col_get(Cat(Isim),'SN').*BB)./(BB + col_get(Cat(Isim),'PEAK_GRADBACK').*InPar.GradBackSize)>InPar.Thresh | PeakForce;
        %FlagBackGrad = (col_get(Cat(Isim),'SN').*BS)./sqrt((BS.^2 + col_get(Cat(Isim),'PEAK_GRADBACK').*InPar.GradBackSize))>InPar.Thresh;
        Cat(Isim)    = row_select(Cat(Isim),FlagBackGrad);
        
        % nead to update Peak list for source flagging
        PeakX = PeakX(FlagBackGrad);
        PeakY = PeakY(FlagBackGrad);
    end    
    
    if (InPar.SearchCR)
        % check if can perform CR
        
        if (~isnan(col_get(Cat(Isim),'PSF_CHI2BACK') ))
        
            for IcrMethod=1:1:numel(InPar.MethodCR)
                % note that the 'chi2backcr' will work only if there are enough
                % stars in the image
                Cat(Isim) = flag_cr_mextractor(Cat(Isim),'Method',InPar.MethodCR{IcrMethod});
            end
        end
    end
   
  
   % Source flags
    if (Calc.Flags)           
        if (InPar.Mask && isfield_populated(Sim(Isim),MaskField))
            % If Mask is true then propogate the Mask into the source Flags
            % find the indices of pixels that belongs to each source
            CellSrcInd = ImUtil.Im.find_within_radius_cell(SizeSimB,PeakX,PeakY,InPar.FlagRadius,true);
            % Calculate the integrated Flag (within FlagRadius) for each
            % source
            Flags = mask4src(Cat(Isim),CellSrcInd,InPar.FlagOperation);
        else
            % do not propogate the Mask
            Nsrc  = size(Cat(Isim).(CatField),1);
            Flags = zeros(Nsrc,1,InPar.FlagClass);
        end
        
        ColFlags = colname2ind(Cat(Isim),'FLAGS');
        Cat(Isim).(CatField)(:,ColFlags) = Flags;
    end
    
    if (InPar.CleanCR)
        % remove CR from source list
        ColFlags  = colname2ind(Cat(Isim),'FLAGS');
        FlagNotCR = ~MASK.isbit_flag(Cat(Isim).(CatField)(:,ColFlags),InPar.CleanBitNames);
        Cat(Isim).(CatField) = Cat(Isim).(CatField)(FlagNotCR,:);
        
    end
    
    if (InPar.CleanDiffSpike)
        % remove sources on diffraction spikes
        FlagNotSpike = ~MASK.isbit_flag(col_get(Cat(Isim),'FLAGS'),InPar.Bit_Spike);
        Cat(Isim).(CatField) = Cat(Isim).(CatField)(FlagNotSpike,:);
    end
    
%     if (InPar.CleanDiffSpike)
%         % remove sources on diffraction spikes
%         FlagSpike = ~MASK.isbit_flag(col_get(Cat(Isim),'FLAGS'),InPar.Bit_Spike);
%         Cat(Isim).(CatField) = Cat(Isim).(CatField)(FlagSpike,:);
%         % nead to update Peak list for source flagging
%         PeakX = PeakX(FlagSpike);
%         PeakY = PeakY(FlagSpike);
%     end  
    
    
    if (~isempty(InPar.SortOut))
        ColIndSort = colname2ind(Cat(Isim),InPar.SortOut);
        if (~isempty(ColIndSort))
            Cat(Isim) = sortrows(Cat(Isim),ColIndSort);
        end  
    end
end
    

% % Set output type
% Cat = colcell2col(Cat);
% switch lower(InPar.OutType)
%     case 'sim'
%         Cat = astcat2sim(Cat,Sim);
%     case 'astcat'
%         % do nothing
%     otherwise
%         error('Unknown OutType option');
% end





end

%--------------------------
%--- Auxilary functions ---
%--------------------------

%---
function Cell=add_index2cell(Cell)
%
Ncell = numel(Cell);
for Icell=1:1:Ncell
    Cell{Icell} = sprintf('%s_%d',Cell{Icell},Icell);
end

end

%---
function Cell=expand_multiple_cell(Cell,Cols,Nrep)
%
Nac       = numel(Cols);
for Iac=1:1:Nac
    Ind = find(strcmp(Cell,Cols{Iac}));
    % duplicate the aperture phot columns
    if (~isempty(Ind))
        Cell = Util.cell.cell_insert(Cell,Cell(Ind),Ind.*ones(1,Nrep-1));
    end
end

end
