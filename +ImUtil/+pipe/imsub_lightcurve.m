function [Res]=imsub_lightcurve(RA,Dec,varargin)
% Generate image subtraction light curve from images
% Package: ImUtil.pipe
% Description: Given images of a transient locatio, read the images,
%              register them, construct a reference image, apply image
%              subtraction (Zackay, Ofek & Gal-Yam 2017) and construct an
%              image subtraction light curve at the transient location.
% Input  : - RA [radians] or NED source name.
%          - Dec [radians].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'NameServerFun' - name server function handle.
%                       Default is @VO.name.server_ned.
%            'DataDir'- directory at which the images resides.
%                       If empty use current directory.
%                       Default is empty.
%            'GetPTF' - Get PTF images for the requested coordinates.
%                       Default is false.
%            'CleanFiles' - Delete images after reduction.
%                       Default is false.
%            'Filter' - A cell array containing a list of filters to
%                       retrieve and analyze.
%                       Default is {'g','R'}.
%            'DirPerFilter' - A flag indicating if there is (or to create)
%                       a directory per filter (true) or to work in a
%                       single directory (false).
%                       Default is true.
%            'JD'     - JD range to retrieve.
%                       Default is [-Inf Inf].
%            'RefJD'  - A matrix of JD range from which to construct a
%                       reference image. The first column is for start JD
%                       and the second column is for end JD.
%            'SubCoaddFlux' - The program calculate the median value of the
%                       flux at the transient location as measured in the
%                       subtraction images from which the reference was
%                       constructed. If true then this flux will be
%                       subtracted from the LC.
%                       Default is true.
%                       This value is stored in the resulted structure
%                       array. Note that this value may be different than
%                       zero due to the fact that the background is
%                       position dependent.
%            'SectionHalfSize' - Half size [X Y] of image cutout on which
%                       the subtraction will be done.
%                       Default is [511 511].
%            'Buffer' - Number of pixels to remove from image edges after
%                       alignment. Default is 3.
%            'MinDist'- The minimum distance of the transient location from
%                       the image edges.  Images in which the transient
%                       distance from image edge is smaller than this value
%                       will not be used. Default is 64.
%            'ExtraHalfSize' - Number of pixels to add to image section
%                       cutouts. Default is 3.
%            'BackPar'- A cell array of additional arguments to pass to the
%                       background.m function.
%                       Default is {'Block',[256 256]}.
%            'MaxBackLevel' - Discard images which mean background level is
%                       higher than this value.
%                       Default is 5000 electrons.
%            'MinNstar' - Discard images in which the number of extracted
%                       sources is smaller than this value.
%                       Default is 20.
%            'MaxAsymptoticRMS' - Discard images in which the astrometric
%                       registration asymptotic RMS is larger than this
%                       value. Default is 0.1 pix.
%            'FlagSatPar' - A cell array of additional arguments to pass to
%                       the flag_saturated SIM method.
%                       Default is {}.
%            'FlagHolesPar' - A cell array of additional arguments to pass
%                       to the flag_holes SIM method.
%                       Default is {}.
%            'AlignPar' - A cell array of additional arguments to pass to
%                       the align SIM method.
%                       Default is {}.
%            'MagCol' - Magnitude column name in the AstCat object.
%                       Default is 'MAG_PSF'.
%            'MagErrCol'- Magnitude error column name in the AstCat object.
%                       Default is 'MAGERR_PSF'.
%            'xcatPar'  - A cell array of additional arguments to pass to
%                       the xcat function.
%                       Default is {}.
%            'RefCat'   - Photometric calibration reference catalog name.
%                       Default is 'sdss'. See 'xcat' for options.
%            'CalibRadius'- Maximum distance of stars from the transient
%                       to use as as photometric calibrators.
%                       Default is 500 arcsec.
%            'ZP'       - Photometric zero point. Default is 22 mag.
%            'MomRadius'- Moments calculation radius.
%                       Default is 4 pix.
%            'MomSigma' - Moment calculation sigma of Gaussian weight.
%                       Default is 1.5 pix.
%            'MomMaxIter' - Maximum number of moments iterations.
%                       Default is 3.
%            'NrandPos'
%            'Verbose' - Verbose. Default is true. 
% Output : - A structure array with the extracted image subtraction light
%            curve and additional information.
%          * Additional image subtraction data and images.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,ResSim]=ImUtil.pipe.imsub_lightcurve
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

ImageField   = SIM.ImageField;
BackField    = SIM.BackField;
ErrField     = SIM.ErrField;
CatField     = AstCat.CatField;
ColCellField = AstCat.ColCellField;

%[RA,Dec]=VO.name.server_ned('SN2010jl','r')

DefV.ImageFileTemplate    = '*_sciimg.fits';
DefV.NameServerFun        = @VO.name.server_ned;
DefV.DataDir              = []; %'/raid/eran/projects/PTF/Science/indiv/PTF10aaxf/LCzogy/';  %[];  %'/raid/eran/tmp/';
DefV.GetPTF               = false;
DefV.CleanFiles           = false;
DefV.Filter               = {'g','r','i'};
DefV.FilterHead           = {'zg','zr','zi'};
DefV.DirPerFilter         = true;
DefV.MethodZP             = 'header';
DefV.KeyZP                = 'ZPAVG';   % ZP header keyword 
DefV.JD                   = [-Inf Inf];
DefV.RefJD                = [-Inf celestial.time.julday([1 7 2010])];
DefV.SubCoaddFlux         = true;
DefV.SatLimit             = 50000;
DefV.SectionHalfSize      = [511 511];
DefV.Buffer               = 9;
DefV.MinDist              = 64;
DefV.ExtraHalfSize        = 3;
DefV.PixScale             = 1;
DefV.ImageSize            = [1000 1000];
DefV.BackPar              = {'Block',[256 256]};
DefV.MaxBackLevel         = 5000;
DefV.MinNstar             = 20;
DefV.MaxAsymptoticRMS     = 0.1;
DefV.FlagSatPar           = {};
DefV.FlagHolesPar         = {};
DefV.AlignPar             = {};
%DefV.MinNumberCommonStars = 10;
DefV.MagCol               = 'MAG_PSF';
DefV.MagErrCol            = 'MAGERR_PSF';
DefV.xcatPar              = {};
DefV.RefCat               = 'SDSSDR10';
DefV.CalibRadius          = 500;   % arcsec
DefV.ZP                   = 22;
DefV.AstromRMS            = 'rrms'; %'asymptotic';  % asymptotic | rrms
DefV.MethodBeta           = 'zp';  % 'zp' | 'zogy'
%DefV.MexThresh            = 7;
DefV.MomRadius            = 4;
DefV.MomSigma             = 1.5;
DefV.MomMaxIter           = 3;
DefV.NrandPos             = 100;
DefV.Verbose              = true;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% If RA is an object name - get RA/Dec from NED name server
if (ischar(RA))
    [RA,Dec] = InPar.NameServerFun(RA,'r');
end

% make sure the filter info is in a cell array
if (~iscell(InPar.Filter))
    InPar.Filter = {InPar.Filter};
end
if (~iscell(InPar.FilterHead))
    InPar.FilterHead = {InPar.FilterHead};
end


% cd to the data directory
PWD = pwd;
if ~isempty(InPar.DataDir)
    cd(InPar.DataDir);
end


Nband = numel(InPar.Filter);
% for each band
for Iband=1:1:Nband
    % Band name is the usual band name (e.g., 'r')
    Band     = InPar.Filter{Iband};
    % BandHead is the band name as appear in the header (e.g., 'zr')
    BandHead = InPar.FilterHead{Iband};
    
    if (InPar.Verbose)
        fprintf('Download and analyze %s-band images\n',Band);
    end
    
    % cd to Band directory
    if (InPar.DirPerFilter)
        if (~exist(BandHead,'dir'))
            mkdir(BandHead);
        end
        cd(BandHead);
    end
    
    % Retrieve images from PTF
    if (InPar.GetPTF)
        [Link,TableSort,Col,FileName]=VO.PTF.wget_corrim(RA,Dec,'Filter',Band,'JD',InPar.JD,'save','i');
        if (isempty(FileName))
            Files = {};
        else
            Files = FileName.Files_Image;
        end
    end
    
    % check for files in dir
    FilesD = dir(InPar.ImageFileTemplate);
    Files = {FilesD.name};
    
    % number of files
    Nim = numel(Files);
    
    
    if (InPar.Verbose)
        fprintf('%d files identified\n',Nim);
    end
    
    if (Nim>0)
        % read header
        % select overlapping X/Y regions that covers RA/Dec position.
        % remove images in which the requested RA/Dec is near boundry
        % Read image sections to SIM object
        % correct for gain
        % Sim local background
        % Sim globacl background
        % remove images with very high background level
        % Extract sources for alignment
        [Sim,SimGlobalBack,Files,Section,CenterSection]=ImUtil.pipe.read_images_section_around_coo(Files,RA,Dec,...
                                                            'MexPar',{},...
                                                            'ZP',InPar.ZP,...
                                                            'MinNstar',InPar.MinNstar,...
                                                            'MaxBackLevel',InPar.MaxBackLevel,...
                                                            'BackPar',InPar.BackPar,...
                                                            'SectionHalfSize',InPar.SectionHalfSize,...
                                                            'ExtraHalfSize',InPar.ExtraHalfSize,...
                                                            'MinDist',InPar.MinDist,...
                                                            'Verbose',InPar.Verbose);
                                                        
        
        % run SWarp without background subtraction
        [FileName,R,AlSim] = ImUtil.pipe.register_swarp(Sim,RA,Dec,'PixScale',InPar.PixScale,'ImageSize',InPar.ImageSize);
        
        % populate the catalog
        AlSim = mextractor(AlSim,'ZP',InPar.ZP);
        
        % select best image
        % image with largest number of stars
        Ncat = sizecat(AlSim);
        [~,IndImMaxStars] = max(Ncat);
        
        
        % JD of images
        JD = julday(Sim);
        AM = mgetkey(Sim,'AIRMASS');
        AM = cell2mat(AM);    

        
        % ZP of images
        switch lower(InPar.MethodZP)
            case 'header'
                % Get ZP from image headers
                ZP = mgetkey(Sim,InPar.KeyZP);
                ZP = cell2mat(ZP);
                
                % flux ZP (transperancy)
                Fzp = 10.^(-0.4.*(26-ZP));
                Fzp = Fzp./mean(Fzp);
                
            case 'calc'
                %--- FIX ---
                
                % Match catalogs and estimate relative ZP (relative to IndRef
                % image)
                [ZP1,Fzp1,ResZP1]=zp_estimate(AlSim,IndImMaxStars,'SkipWCS',true,...
                                                      'Columns',{'XWIN_IMAGE','YWIN_IMAGE',InPar.MagCol,InPar.MagErrCol},...
                                                      'MagCol',InPar.MagCol,...
                                                      'MagErrCol',InPar.MagErrCol);
        
            otherwise
                error('Unknown MethodZP option');
        end
        
        
        %--- Prep images to coaddition ---
        % subtract local background, but keep global background!
        % keep global background in the background field
        % note there is a problem near the bright galaxy - likely related to background subtraction...
        % therefor use global background

        Nsim = numel(AlSim);
        for Isim=1:1:Nsim
            AlSim(Isim) = sub_background(AlSim(Isim));
            AlSim(Isim).Im = AlSim(Isim).Im + SimGlobalBack(Isim).BackIm;
            AlSim(Isim).BackIm = SimGlobalBack(Isim).BackIm;
            AlSim(Isim).ErrIm  = SimGlobalBack(Isim).ErrIm;
        end

        % repopulate Mask:
        % in the future - change scheme
        AlSim = flag_holes(AlSim,'Nsigma',3,InPar.FlagHolesPar{:});
        AlSim = flag_saturated(AlSim,InPar.FlagSatPar{:});
        AlSim = flag_cr_mextractor(AlSim);
        AlSim = flag_cr_mextractor(AlSim,'Method','chi2backcr');
        %AlSim = flag_backgrad(AlSim);
        Nim = numel(AlSim);
        
        
        % select images for reference image:
        SelIndRef = Util.array.find_ranges(JD,InPar.RefJD);

        if (InPar.Verbose)
            fprintf('Construct reference image from %d images\n',numel(SelIndRef));
        end

        if (numel(SelIndRef)==0)
            error('Cant find images for construction of Reference');
        end
        
        %--- Coadd images ---
        
        % subtract background and scale images
        % this is needed only for the reference image
        AlSimBS = sub_background(AlSim(SelIndRef));
        % scale the image to the mean level of the images in the ref image
        AlSimBS = scale(AlSimBS,Fzp(SelIndRef));
        NimRef  = numel(AlSimBS);

        % regular coaddition
        [CoaddSim,StdSim,WeightSim]=coadd(AlSimBS,'CoaddMethod','meanclip','MethodZero','none','MethodScale','none');
        % sum the images in order to get an estimate for the background
        [SumCoaddSim]=coadd(AlSim(SelIndRef),'CoaddMethod','sum','MethodZero','none','MethodScale','none');
        SumCoaddSim = background(SumCoaddSim,'Block','full');
        
        % Set the flux level of CoaddSim to be roughly the photon counts
        CoaddSim.(ImageField) = CoaddSim.(ImageField).*NimRef + SumCoaddSim.(BackField);
        
        CoaddSatLimit = InPar.SatLimit.*NimRef + SumCoaddSim.(BackField);
        
        PSFSelectorPar       = {'ColSN','SN',...
                           'ColMaxFlux','PEAKF_VALTOT',...
                           'ColBitFlag','FLAGS',...
                           'MaskVal',0,...
                           'PosCol',{'XWIN_IMAGE','YWIN_IMAGE'},...
                           'BoundryDist',10,...
                           'MinSN',10,...
                           'SatLevel',CoaddSatLimit};
        % find sources in CoaddSim
        CoaddSim = mextractor(CoaddSim,'ZP',InPar.ZP,'PSFSelectorPar',PSFSelectorPar);
        
        % set background to global value
        CoaddSim.(BackField)  = SumCoaddSim.(BackField);
        CoaddSim.(ErrField)   = SumCoaddSim.(ErrField);
        
        
        %--- add astrometric solution to reference image ---
        [CoaddSim_ResAst,CoaddSim] = astrometry(CoaddSim,'RA',RA,'Dec',Dec,'Flip',[1 1;1 -1;-1 1;-1 -1]);   % No flip as SWarp fliped the images!!
        % update RA/Dec in CoaddSim
        CoaddSim = update_coordinates(CoaddSim);
        
        %--- find photometric solution of CoaddSim (Reference image) ---
        ResPhot = ImUtil.Im.photometric_calibration(CoaddSim);
        %ResPhot = ImUtil.Im.photometric_calibration(Sim(IndImMaxStars));

        
        % Target coordinates in the CoaddSim image
        % assume target is in the image center!
        W = ClassWCS.populate(CoaddSim);
        %[TargetX, TargetY] = coo2xy(W,[RA  Dec]);
        TargetX = InPar.ImageSize(1).*0.5;
        TargetY = InPar.ImageSize(2).*0.5;
        
        
        %%
        % estimate the astrometric residual relative to the reference image
        % can be made faster by just calculating residuals 

        % Match catalogs and estimate ZP relative to CoaddSim
        IndRefC = numel(AlSim) + 1;
        [ZPC,FzpC,ResZPC]=zp_estimate([AlSim(:);CoaddSim],IndRefC,'SkipWCS',true,...
                                              'Columns',{'XWIN_IMAGE','YWIN_IMAGE',InPar.MagCol,InPar.MagErrCol},...
                                              'MagCol',InPar.MagCol,...
                                              'MagErrCol',InPar.MagErrCol);
        

        % CoaddSim ZP     
        CoaddImageZP = InPar.ZP - ResPhot.Par(1); % -2.5.*log10(FzpC(IndRef)) ;
        
        % fix PSF for image subtraction
        AlSim    = psf_extrapolate(AlSim);
        CoaddSim = psf_extrapolate(CoaddSim);

        
        switch lower(InPar.AstromRMS)
            case 'asymptotic'
                SigmaX = sqrt([R.MinAssymErr].'.^2 + CoaddSim_ResAst.MinAssymErr.^2) .* 3600./InPar.PixScale;
                SigmaY = SigmaX;
                
                FlagSigN = isnan(SigmaX) | isnan(SigmaY);
                if (any(FlagSigN))
                    SigmaXr = sqrt([R.rrms1].'.^2 + CoaddSim_ResAst.rrms1.^2).*3600./InPar.PixScale;
                    SigmaYr = SigmaXr;
                    SigmaX(FlagSigN) = SigmaXr(FlagSigN);
                    SigmaY(FlagSigN) = SigmaYr(FlagSigN);
                end
            case 'rrms'
                SigmaX = sqrt([R.rrms1].'.^2 + CoaddSim_ResAst.rrms1.^2).*3600./InPar.PixScale;
                SigmaY = SigmaX;
            otherwise
                error('Unknown AstromRMS option');
        end
        
        %AlSim = single(AlSim);
        switch lower(InPar.MethodBeta)
            case 'zp'
%                 [Summary,D,S,Scorr,SigmaF]=subtract_proper(AlSim,CoaddSim,'BackReCalc',false,'FluxMatch',FzpC(1:end-1),...
%                                     'SigmaX',SigmaX,'SigmaY',SigmaY);
                                
                                
                 FluxMatch = Fzp.*mean(FzpC(1:end-1))./mean(Fzp);
                 [Summary,D,S,Scorr,SigmaF]=subtract_proper(AlSim,CoaddSim,'BackReCalc',false,'FluxMatch',FluxMatch,...
                                    'SigmaX',SigmaX,'SigmaY',SigmaY);
                   
                                

                
%                 B1=timeseries.binning([JD-2458000.5, Mag1],1,[150 290],{'MidBin',@mean,@std,@numel});
                 

            case 'zogy'
                [Summary,D,S,Scorr,SigmaF]=subtract_proper(AlSim,CoaddSim,'BackReCalc',false,...
                                    'SigmaX',SigmaX,'SigmaY',SigmaY);
            otherwise
                error('Unknown MethodBeta option');
        end
        
        
        S = background(S);
        
        % Generate LC from subtraction images
        [DataT,DataR] = ImUtil.pipe.lightcurve_from_sub_im([TargetX, TargetY],S,D,Scorr,SigmaF,Summary,CoaddSim,...
                                            'ZP',CoaddImageZP,'SigmaX',SigmaX,'SigmaY',SigmaY);
        
        
        
 
        % FluxRand statistics
        %FluxStdRand = Util.stat.rstd(FluxRand,2);
        % FluxStdRand = std(FluxRand,[],2);

        % Mask
        Filter=Kernel2.aper(7);
        Filter(Filter>0) = 1;
        Filter = uint32(Filter);

        AlSimMask = ufun2sim(AlSim,@imdilate,'ExecField',{'Mask'},'FunAddPar',{Filter});
        AlSimMask = mask_add(AlSimMask,D);

        [~,~,ValMask]=get_value(AlSimMask,[TargetX, TargetY]);
    

        % mean flux offset of the S image
        OffsetS = ([S.BackIm]./[Summary.F_S]).';
        Beta    = [Summary.Beta].';
        GeodPos = geodpos(Sim);
        [AM,AzAlt] = celestial.coo.airmass(JD,RA,Dec,[[GeodPos.Long].',[GeodPos.Lat].']);

        Res(Iband).ZP = CoaddImageZP;
        Res(Iband).RefJD = InPar.RefJD;
        Res(Iband).ColorTerm = ResPhot.Par(2);
        %Res(Iband).RefCatBand = CatBand;
        %Res(Iband).RefCatBandColor = CatColor;
        Res(Iband).Filter = Band;
        Res(Iband).DataT  = DataT;
        Res(Iband).DataR  = DataR;
        Res(Iband).CatLC = AstCat;
        ONES = ones(size(JD));
        Res(Iband).CatLC.(CatField) = [JD, DataT.MagBest, DataT.MagErrBest, DataT.MagErrAllBest, ...
                                           DataT.ValBest, DataT.ValErrBest, DataT.ValErrAllBest, ...
                                           ValMask, ...
                                           DataT.BestX,   DataT.BestY, ...
                                           DataT.ValTar,  DataT.ValTarSc, ...
                                           DataT.RandNon0_Nsig, ...
                                           DataT.RatioMeanF_S, ...
                                           DataR.Nrand.*ONES, DataR.Nepoch.*ONES, ...
                                           DataR.EpochChi2, ...
                                           DataT.F_S, OffsetS, Beta, AM, AzAlt.*RAD];
        Res(Iband).CatLC.(ColCellField) = {'JD', 'Mag','MagErr', 'MagErrAll', ...
                                                 'Flux', 'FluxErr', 'FluxErrAll', ...
                                                 'Mask',...
                                                 'BestX', 'BestY', ...
                                                 'FluxTarget', 'ScTarget', ...
                                                 'RandNon0_Nsig', 'RatioMeanF_S' ...
                                                 'Nrand', 'Nepoch', ...
                                                 'EpochChi2',...
                                                 'F_S','OffsetS','Beta','AM', 'AZ', 'Alt'};
            
        %FluxMM, SigmaValPeakMM, FluxStdRand, SigmaAll, Mag, OffsetS, Beta, AM, AzAlt.*RAD, MomData, FluxCoaddMedian.*ONES, FluxCoaddStd.*ONES, FluxCoaddNim.*ONES, SigmaX(:), SigmaY(:), SGPeak, RandChi2, RandDof.*ONES, ValMask];
        %Res(Iband).CatLC.(ColCellField) = {'JD','FluxMM','FluxErr','FluxStdRand','SigmaAll','Mag','OffsetS','Beta','AM','Az','Alt','Xforce','Yforce','X','Y','X2','Y2','XY','Xmax','Ymax','FluxCoaddMedian','FluxCoaddStd','FluxCoaddNim','SigmaX','SigmaY','GradPeak','RandChi2','RandDof','FLAGS'};

        Res(Iband).CatLC = colcell2col(Res(Iband).CatLC);
        Res(Iband).CatLC = sortrows(Res(Iband).CatLC,1);

        
        

%         % find transients
%         Thresh = 5;
%         Nsim = numel(Scorr);
%         for Isim=1:1:Nsim
%             Icand = find(abs(Scorr(Isim).(ImageField))>Thresh);
%             [Yc,Xc] = ind2sub(size(Scorr(Isim).(ImageField)),Icand);
%             [~,~,ValMask]=get_value(AlSimMask(Isim),[Xc,Yc]);
%             ValMask = squeeze(ValMask);
%             Im0 = find(ValMask==0);
%             Xc = Xc(Im0);
%             Yc = Yc(Im0);
%             [TMom(Id),TMom2(Id)]=ImUtil.Im.im_moments(D(Isim).(ImageField),Xc, Yc,InPar.MomRadius,InPar.MomSigma,InPar.MomMaxIter);
%             [Xc,Yc] %, TMom2(Id).X2, TMom2(Id).Y2, TMom2(Id).XY ]
% 
%             ds9(AlSim(Isim),1); pause(0.5);
%             ds9(Scorr(Isim),2); pause(0.5);
%             ds9(double(abs(Scorr(Isim))>5),3); pause(0.5);
%             ds9(D(Isim),4); pause(0.5);
%             %ds9.lock_xy(true); pause(0.5);
%             ds9.plot(Xc,Yc);
%         end



        if (InPar.CleanFiles)
            % delete PTF images
            Util.files.delete_cell(FileName);
        end
    end
    
    if (InPar.DirPerFilter)
        cd('../');
    end
    
%     % write results to file:
    TextFile = sprintf('LC_Band_%s.txt',Res(Iband).Filter);
    FID      = fopen(TextFile,'w');
    fprintf(FID,'%% Generated by: ImUtil.pipe.imsub_lightcurve\n');
    fprintf(FID,'%% https://webhome.weizmann.ac.il/home/eofek/matlab/index.html\n');
    fprintf(FID,'%% Written by Eran Ofek\n');
    fprintf(FID,'%% \n');
    fprintf(FID,'%% If you use this data in scientific work - please provide the following credit:\n');
    fprintf(FID,'%%     The light curve was generated by an image subtraction light curve package\n');
    fprintf(FID,'%%     available as part of the MATLAB Astronomy and Astropgyscis Toolbox (Ofek 2014),\n');
    fprintf(FID,'%%     using the  proper image subtraction described in Zackay, Ofek, Gal-Yam (2016).\n');
    fprintf(FID,'%%     The reference image was constructed using the weights discussed in.\n');
    fprintf(FID,'%%     Zackay and Ofek (2017a, 2017b). The images astrometry was conducted using\n');
    fprintf(FID,'%%     {\tt mextractor} (Ofek 2018), and {\tt astrometry.m} (Ofek 2018),\n');
    fprintf(FID,'%%     while the image registration was done using {\tt SWarp} (Bertin 2010).\n');
    fprintf(FID,'%% \n');
    fprintf(FID,'%% Generation date: %04d-%02d-%02d %02d:%02d:%04.1f\n',clock);
    fprintf(FID,'%% \n');
    fprintf(FID,'%% Object J2000.0 R.A.: %11.6f\n',RA.*RAD);
    fprintf(FID,'%% Object J2000.0 Dec.: %11.6f\n',Dec.*RAD);
    fprintf(FID,'%% Filter: %s\n',Band);
    %fprintf(FID,'%% Number of images in reference: %d\n',Res(Iband).FluxCoaddNim);
    fprintf(FID,'%% Reference image JD: %12.2f -- %12.2f\n',Res(Iband).RefJD);
    fprintf(FID,'%% ZP: %8.3f\n',Res(Iband).ZP);
    fprintf(FID,'%% ColorTerm: %8.3f\n',Res(Iband).ColorTerm);
    fprintf(FID,'%% ');
    for Icol=1:1:numel(Res(Iband).CatLC.ColCell)
        fprintf(FID,'%s  ',Res(Iband).CatLC.ColCell{Icol});
    end
    fprintf(FID,'\n');
    fprintf(FID,'%14.6f %8.3f %8.3f %8.3f  %+10.5e %+10.5e %+10.5e  %6d %8.2f %8.2f %+10.5e %7.2f %7.2f %7.3f %6d %6d %10.2f %+10.5e %+10.5e %+10.5e %6.3f %6.2f %6.2f\n',Res(Iband).CatLC.Cat.');
    fclose(FID);
    
    
end

% cd to original directory
cd(PWD);