function [OutPSF,Res]=psf_estimator(Sim,varargin)
% PSF estimator.OBSOLETE: Use psf_extractor instead.
% Package: @SIM
% Description: Estimate the PSF of a SIM object.
%              OBSOLETE: Use psf_extractor instead.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ImageField' - The image field name in the SIM, from which to
%                           estimate the PSF. Default is 'Im'.
%            'BackField'  - The background image field in the SIM.
%                           Default is 'BackIm'.
%            'ErrField'   - The error image field in the SIM.
%                           Default is 'ErrIm'.
%            'Align'      - Align the PSF stamps. Default is true.
%                           See stamp_xy.m for details.
%            'HalfPixShift'-Shift the stamp by additional half pixel.
%                           This controls if the source center is in the
%                           pixel center or on a pixel edge.
%                           Default is false.
%                           The default option (false) will locate the
%                           center of the PSF in the middle of the central
%                           pixel, where the central pixel is given by
%                           ceil(N./2+0.1).
%            'Gain'       - Gain [e-/ADU] of images.
%                           If numeric, then apply as Gain.
%                           If string or cell array of strings than use
%                           getkey_gain.m to get the image Gain from the
%                           header.
%                           If empty, then use getkey_gain.m to get the
%                           image Gain from the header using its default
%                           parameters. Default is [].
%            'SubBack'    - Subtract background from images.
%                           Default is true.
%            'BackSim'    - A SIM array of background images in the 'Im'
%                           field (in e- units!). If not provided then
%                           attempt to use the take the background from
%                           the 'BackIm' field (Gain will be applied).
%                           Default is [].
%            'ReEstBack'  - Re-estimate background even if provided in
%                           the 'BackIm' field. Default is false.
%            'BackPar'    - Cell array of arguments to pass to the
%                           background.m function.
%                           Default is
%                           {'MethodBack','mode_fit','MethodStD','mode_fit','Block',[256 256]}.
%            'PosCat'     - An AstCat object, or a two column matrix,
%                           containing the X,Y positions of
%                           the sources in the image. If this is a matrix
%                           then assume only two columns are provided
%                           [X,Y].
%                           If empty then will attempt to read the catalog
%                           from the SIM object. If empty and not SIM
%                           catalog entry then will attempt to repopulate
%                           the catalog (See 'FunCat' parameter).
%                           Default is empty.
%            'PosCatCol'  - Cell array of columns containing the X and Y
%                           position columns in the 'PosCat' catalog.
%                           Default is {'XWIN_IMAGE','YWIN_IMAGE'}.
%            'MomCatCol'  - Cell array of columns containing the X^2, Y^2
%                           and X*Y moment columns in the 'PosCat' catalog.
%                           Default is
%                           {'X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE'}.
%            'ReCat'      - If 'PosCat' is empty, then this is a logical
%                           flag indicating if to re-populate the
%                           'PosCat' catalog (in the SIM input) using the
%                           function specified in 'FunCat'.
%                           Default is false.
%            'FunCat'     - Function to use for poulating the source
%                           catalog. Default is @sextractor.
%            'FunCatPar'  - Cell array of arguments to pass to 'FunCat'.
%                           Default is {}.
%            'SelectFun'  - PSF candidate selection function.
%                           A function that select the best candidate PSF
%                           stars in a acatalog. The function input and
%                           output are Astcat objects.
%                           AstCat=Fun(AstCat,SelectFunPar{:}).
%                           See an example in psf_cat_selector.m.
%                           Default is @psf_cat_selector.
%            'SelectFunPar'- A cell array of arguments to pass to
%                           'SelectFun'. Default is {}.
%            'StampSize'   - PSF stamp half size. Default is 10.
%            'MomentSigma' - If the second moments are not provided in the
%                            catalog input then this is the Gaussian sigma
%                            used for the window weighting in the moments
%                            calculation. Default is 1.5 pix.
%            'UseErr'      - A logical flag indicating if to use for the
%                            background variance the SIM image 'ErrIm'
%                            field squared (true), or the image local
%                            background (false). Default is true.
%            'MethodChi2Thresh'- A mrthod by which to select PSF stars
%                            based on their \chi^2 relative to the mean
%                            PSF. Options are:
%                            'prob' - Use chi2inv.m threshold.
%                                     In this case 'ParChi2Thresh' contains
%                                     the probability (e.g., 0.999).
%                            'min'  - Multiple the minimum \chi^2 by a
%                                     factor given by 'ParChi2Thresh'
%                                     (e.g., 10).
%                            'quant'- Choose the 'ParChi2Thresh' lower
%                                     quantile of best \chi^2 values
%                                     (e.g., 0.3).
%                            Default is 'quiant'.
%            'ParChi2Thresh'- Parameters for 'MethodChi2Thresh'.
%                            Default is 0.3.
%            'NoiseFactor'  - The minimum S/N to accept as valid measurment
%                            for the mean PSF. All the PSF pixels, with 
%                            radii larger than the minimum radius of pixel
%                            which S/N is lower than this value will be
%                            set to 0. Default is 8.
%            'Rad0'        - Set the PSF values above this radius to 0.
%                            Default is 6.
%            'NoiseEstimate'- Method by which to estimate the PSF error.
%                            Options are:
%                            'grad' - Using the PSF gradient.
%                            'wmean'- Using weighted mean. Default.
%            'OutType'      - One of the following output types:
%                             'OrigSim' - PSF object is added to the
%                                         original input SIM object.
%                             'NewSim'  - PSF object is added to the
%                                         modified input SIM object.
%                             'PSF'     - Only a PSF object.
%                             Default is 'origsim'.
%            'Verbose'      - Verbose. Default is false.
% Output : - The estimated PSF.
%            See 'OutType' for output options.
%          - Structure array of statistics regarding the estimated PSF.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [OutPSF,Res]=psf_estimator(S1);
% Reliable: 2
%--------------------------------------------------------------------------
import timeseries.*

ImageField = 'Im';
BackField  = 'BackIm';
ErrField   = 'ErrIm';
CatField   = 'Cat';


DefV.ImageField         = ImageField;
DefV.BackField          = BackField;
DefV.ErrField           = ErrField;
DefV.Align              = true;
DefV.HalfPixShift       = false;
DefV.Gain               = [];   % [] read from header; Number; Key name
DefV.SubBack            = true;
DefV.BackSim            = [];   % stored in the ImageField field
DefV.ReEstBack          = false;
DefV.BackPar            = {'MethodBack','mode_fit','MethodStD','mode_fit','Block',[256 256]};
DefV.PosCat             = [];   % matrix or AstCat of PSF stars positions
DefV.PosCatCol          = {'XWIN_IMAGE','YWIN_IMAGE'};
DefV.MomCatCol          = {'X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE'};
DefV.ReCat              = false;
DefV.FunCat             = @mextractor;
DefV.FunCatPar          = {};
DefV.SelectFun          = @psf_cat_selector;
DefV.SelectFunPar       = {};
DefV.StampSize          = 10;
DefV.MomentSigma        = 1.5;   % sigma used for weighting the moments calculation
DefV.UseErr             = true;
DefV.MethodChi2Thresh   = 'quant';
DefV.ParChi2Thresh      = 0.7;
DefV.NoiseFactor        = 3;
DefV.Rad0               = 6;
DefV.NoiseEstimate      = 'grad'; %'wmean';
DefV.OutType            = 'OrigSim';   % 'NewSim' | 'OrigSim' | 'PSF' |
DefV.Verbose            = false;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Nsim  = numel(Sim);
Nback = numel(InPar.BackSim);

% Allocate output
switch lower(InPar.OutType)
    case 'origsim'
        OutPSF = Sim;
    case 'newsim'
        % do nothing
    case 'psf'
        OutPSF = psfdef(size(Sim));
    otherwise
        error('Unknown OutType option');
end


% Get gain for all images.
if (isempty(InPar.Gain))
    % read Gain from header as default
    Gain = getkey_gain(Sim);
elseif (ischar(InPar.Gain) || iscellstr(InPar.Gain))
    Gain = getkey_gain(Sim,'GainKeys',InPar.Gain);
elseif (isnumeric(InPar.Gain))
    Gain = InPar.Gain.*ones(Nsim,1);
else
    error('Unknown Gain option');
end

if (InPar.Verbose)
    fprintf('Gain range is %f to %f\n',min(Gain),max(Gain));
end

%------------------------
%--- Correct for Gain ---
%------------------------
% Apply Gain correction to image, background and Err
if (InPar.Verbose)
    fprintf('Correct for Gain\n');
end
if (all(Gain==1))
    % do nothing - no Gain correction
else
    Sim = bfun2sim(Sim,Gain,@times,'ExecField',{InPar.ImageField, InPar.BackField, InPar.ErrField});
end


%---------------------------
%--- subtract background ---
%---------------------------
if (InPar.SubBack)
    if (InPar.ReEstBack)
        % re-estimate background - even if provided
        Sim = background(Sim,'ExecField',InPar.ImageField,...
                             'SubBack',InPar.SubBack,...
                             InPar.BackPar{:});
                 
    else
        % background is provided by user, but still not subtracted
        if (isempty(InPar.BackSim))
            % check if background is in Sim(Isim).BackIm:
            FlagBack = Util.struct.isfield_notempty(Sim,BackField);   % BackIm field exist and not empty
            % Estimare background for missing background cases
            % Note that here the background is not subtracted...
            Sim(~FlagBack) = background(Sim(~FlagBack),'ExecField',InPar.ImageField,...
                                            'SubBack',false,...
                                            InPar.BackPar{:});
            
            % assume background is stored in Sim.BackIm
            % subtract background from each image
            for Isim=1:1:Nsim
                Sim(Isim).(InPar.ImageField) = Sim(Isim).(InPar.ImageField) - Sim(Isim).(InPar.BackField);
            end

        else
            % assume background is stored in InPar.BackSim
            for Isim=1:1:Nsim
                Iback = min(Isim,Nback);
                Sim(Isim).(InPar.ImageField) = Sim(Isim).(InPar.ImageField) - InPar.BackSim(Iback).(InPar.ImageField);
            end
        end
    end
end
if (InPar.Verbose)
    fprintf('Background subtracted\n');
end
% Sim is now background subtracted


%-------------------------------------------
%--- Generate/get positions of PSF stars ---
%-------------------------------------------
if (isempty(InPar.PosCat))
    % The user didn't supply the positions for the PSF stars
    
    % generate the associated source catalog
    if (InPar.ReCat)
        % Re-generate source catalog even if exist
        Sim = InPar.FunCat(Sim,InPar.FunCatPar{:});
    else
        % Use existing source catalog if exist
        FlagPopCat = iscat_populated(Sim);
        if ~all(FlagPopCat)
            Sim(~FlagPopCat) = InPar.FunCat(Sim(~FlagPopCat),InPar.FunCatPar{:});   
        end
    end

    % Copy AstCat into its own structure...
    %PosCat = sim2astcat(Sim);
    
    % Use the PosCat in the Sim array
    PosCat = Sim;
    
else
    % The user supplied the final positions for the PSF stars
    if (~AstCat.isastcat(InPar.PosCat))
        % Assume InPar.PosCat is in a two column matrix form
        % convert to AstCat
        % Note that this produce single element catalog that will be use
        % for all SIM elements.
        PosCat = AstCat.array2astcat(InPar.PosCat,InPar.PosCatCol);
    else
        % InPar.PosCat is already in AstCat format - do nothing
        PosCat = InPar.PosCat;
    end
end
if (InPar.Verbose)
    fprintf('Source positions obtained\n');
end

% Allocate output
switch lower(InPar.OutType)
    case 'newsim'
        OutPSF = Sim;
    otherwise
        % do nothing
end


%-----------------------------------------------
%--- Preliminary selection of best PSF stars ---
%-----------------------------------------------
% as minimum make sure that all selected stars are within the image
CCDSEC  = ccdsec(Sim);
CCDSECb = [CCDSEC(:,1) + InPar.StampSize + 1,...
           CCDSEC(:,2) - InPar.StampSize - 1,...
           CCDSEC(:,3) + InPar.StampSize + 1,...
           CCDSEC(:,4) - InPar.StampSize - 1];
[PosCat,~]=select_ccdsec(PosCat,CCDSECb,InPar.PosCatCol);
% select PSF stars
if (~isempty(InPar.SelectFun))
    [PosCat] = InPar.SelectFun(PosCat,InPar.SelectFunPar{:});
else
    % Use the PosCat in the Sim array
    %PosCat = PosCat;
end
if (InPar.Verbose)
    fprintf('Candidate PSF star positions obtained\n');
end

% Define X,Y, radius coordinates for stamp
[MatX,MatY] = meshgrid((-InPar.StampSize:1:InPar.StampSize),(-InPar.StampSize:1:InPar.StampSize));
MatR = sqrt(MatX.^2 + MatY.^2);

%-------------------------------------------------------------
%--- get stamps for all positions - this is done per image ---
%-------------------------------------------------------------
Res  = Util.struct.struct_def({'Cube','Npsf','X','Y','Back','Flux','XstM','YstM','X2stM','Y2stM','XYstM','M2','Chi2','Flag','Ngood','R0','PSF','ErrPSF'},size(Sim));
Npos = numel(PosCat);
for Isim=1:1:Nsim
    if (InPar.Verbose)
        fprintf('Get stamps for image %d\n',Isim);
    end
    % for each SIM image
    Ipos   = min(Isim,Npos);
    
    % construct a cube of all PSF stars aligned
    StampWidth = InPar.StampSize.*2 + 1;
    % Note that the output is cube [star,x,y]:
    if (sizecat(PosCat(Ipos))==0)
        error('No candidate PSF stars were found');
    end
    [Out] = stamp_xy(Sim(Isim),PosCat(Ipos),'ExecField',      InPar.ImageField,...
                                            'StampSize',      InPar.StampSize,...
                                            'Align',          InPar.Align,...
                                            'HalfPixShift',   InPar.HalfPixShift,...
                                            'OutType',        'cube',...
                                            'PadValue',       NaN,...
                                            'PosCatCol',      InPar.PosCatCol);
    
    %----------------------------------------------
    %--- Calculate the properties of the stamps ---
    %----------------------------------------------
    if (InPar.Verbose)
        fprintf('Calculate stamps properties\n');
    end
    Res(Isim).Npsf = size(Out,1);
    % read background at positions
    ColInd = colname2ind(PosCat(Ipos),InPar.PosCatCol);
    XY  = PosCat(Ipos).(CatField)(:,ColInd);
    
    
    Ind = sub2ind(size(Sim(Isim).(InPar.ImageField)),round(XY(:,2)),round(XY(:,1)));
    if (numel(Sim(Isim).(BackField))==1)
        IndB = 1;
    else
        IndB = Ind;
    end
    if (numel(Sim(Isim).(ErrField))==1)
        IndE = 1;
    else
        IndE = Ind;
    end
    
    BackLevel = Sim(Isim).(BackField)(IndB);
    if (InPar.UseErr)
        % Use the local estimated std^2 as the local background variance
        BackVariance  = (Sim(Isim).(ErrField)(IndE)).^2;
    else
        % Use the local estimated background (e- units) as the local
        % background variance.
        % Note that this ignores RN^2 bias and flat errors...
        BackVariance = BackLevel;
    end
    
    % calculate the annulus background for each stamp
    Res(Isim).X    = col_get(PosCat(Ipos),InPar.PosCatCol{1});
    Res(Isim).Y    = col_get(PosCat(Ipos),InPar.PosCatCol{2});
    
    [AperPhot,AperPhotCol] = aper_phot(Sim(Isim),[Res(Isim).X, Res(Isim).Y]);
    BackLevel = BackLevel + AperPhot.(CatField)(:,AperPhotCol.Back);
    Out = bsxfun(@minus,Out,AperPhot.(CatField)(:,AperPhotCol.Back));
    
    
    % Note that in this stage the PSF is assumed to be background subtracted
    Res(Isim).Cube = Out;
    %Res(Isim).X    = col_get(PosCat(Ipos),InPar.PosCatCol{1});
    %Res(Isim).Y    = col_get(PosCat(Ipos),InPar.PosCatCol{2});
    Res(Isim).Back = BackLevel;
    Res(Isim).Flux = sum(sum(Out,2),3);
    Nstar = numel(Res(Isim).X);
    BackVariance = BackVariance.*ones(Nstar,1);
    
    StampCenInd = StampWidth.*0.5 + 0.5;
    
    
    
    % If moments exist in catalog then use them
    % if not then calculate moments
    
    InPar.Select2nMom = true;
    if (InPar.Select2nMom)
        % Select stamps based on second moment distribution as a function
        % of flux
        
        %VecStamp    = (-InPar.StampSize:1:InPar.StampSize);
        if (isastcat_col(Sim(Isim),InPar.MomCatCol{1}) && isastcat_col(Sim(Isim),InPar.MomCatCol{2}) && isastcat_col(Sim(Isim),InPar.MomCatCol{3}) )
            % moment columns exist in catalog
            Res(Isim).X2stM = col_get(PosCat(Ipos),InPar.MomCatCol{1});
            Res(Isim).Y2stM = col_get(PosCat(Ipos),InPar.MomCatCol{2});
            Res(Isim).XYstM = col_get(PosCat(Ipos),InPar.MomCatCol{3});
        else
            % moment columns not in catalog - calculate
            for I=1:1:Res(Isim).Npsf
                %R(I) = moment2d(VecStamp,VecStamp,squeeze(Out(I,:,:)),StampCenInd,StampCenInd);    
                [Mom(I),Mom2(I),~]=ImUtil.Im.im_moments(squeeze(Out(I,:,:)),StampCenInd,StampCenInd,InPar.StampSize,InPar.MomentSigma);
            end
            Res(Isim).XstM  = [Mom.X].';
            Res(Isim).YstM  = [Mom.Y].';
            Res(Isim).X2stM = [Mom2.X2].';
            Res(Isim).Y2stM = [Mom2.Y2].';
            Res(Isim).XYstM = [Mom2.XY].';
        end
        Res(Isim).M2    = Res(Isim).X2stM + Res(Isim).Y2stM;

        %--------------------------
        %--- Select good stamps ---
        %--------------------------
        % select stamp in flux range according to stability of the second
        % moment.
        % I.e., select flux range in which the second moment is stable.
        % This will tend to remove faint and bright stars for which the second
        % moment have large fluctuations from the median value.

        if (InPar.Verbose)
            fprintf('Select good stamps\n');
        end
        %B=binning([log10(Res(Isim).Flux),Res(Isim).M2],0.1);
        B = binning([log10(Res(Isim).Flux),Res(Isim).M2],0.1,[NaN NaN],{'MidBin',@median,@Util.stat.mean_error});
        Bin.LogF = B(:,1); %B(:,1);
        Bin.MedY = B(:,2); %B(:,8);
        Bin.ErrY = B(:,3); %B(:,3);
        NormDiff = diff(Bin.MedY)./sqrt(Bin.ErrY(1:end-1).^2 + Bin.ErrY(2:end).^2);
        %MeanLogF = (Bin.LogF(1:end-1) + Bin.LogF(2:end)).*0.5;
        IndGoodLogF = find(abs(diff(abs(NormDiff)<1)));
        GoodFluxRange = 10.^Bin.LogF(IndGoodLogF(1:2));
        FlagGoodFlux = Res(Isim).Flux>GoodFluxRange(1) & Res(Isim).Flux<GoodFluxRange(2);
        % Flag contains the logical indicating if the star is a good PSF
        % candidate
        Res(Isim).Flag = FlagGoodFlux;
    else
        Res(Isim).Flag = true(size(Res(Isim).X));
    end
    
    %------------------------------
    %--- Construct the mean PSF ---
    %------------------------------
    
    P=squeeze(sum(Res(Isim).Cube(Res(Isim).Flag,:,:),1));
    P = P./sum(P(:));
    ScaledP = bsxfun(@times,Res(Isim).Flux,reshape(P,1,StampWidth,StampWidth));
    Resid = bsxfun(@minus,Res(Isim).Cube,ScaledP);
    % Note the use of BackVariance
    % BackVariance is either background level, or background std^2
    % (according to UseErr).
    Chi2 = Resid.^2./bsxfun(@plus,Res(Isim).Cube,BackVariance);
    Chi2 = sum(sum(Chi2,2),3);
    
    
    % Set the Chi2 threshold for selection of good PSF candidates
    switch lower(InPar.MethodChi2Thresh)
        case 'prob'
            ThreshChi2 = chi2inv(InPar.ParChi2Thresh,StampWidth.^2);
        case 'min'
            % select ParChi2Thresh times the min(Chi2)
            ThreshChi2 = min(Chi2).*InPar.ParChi2Thresh;
        case 'quant'
            % select ParChi2Thresh lower quantile
            ThreshChi2 = quantile(Chi2,InPar.ParChi2Thresh);
        otherwise
            error('Unknown MethodChi2Thresh option');
    end
    Res(Isim).Flag = Chi2<ThreshChi2 & Res(Isim).Flag;
    
    % second iteration PSF
    P=squeeze(sum(Res(Isim).Cube(Res(Isim).Flag,:,:),1));
    P = P./sum(P(:));
    ScaledP = bsxfun(@times,Res(Isim).Flux,reshape(P,1,StampWidth,StampWidth));
    Resid = bsxfun(@minus,Res(Isim).Cube,ScaledP);
    Chi2 = Resid.^2./bsxfun(@plus,Res(Isim).Cube,BackVariance);
    Chi2 = sum(sum(Chi2,2),3);
    
    Res(Isim).Chi2 = Chi2;
    Res(Isim).Ngood = sum(Res(Isim).Flag);
    switch lower(InPar.NoiseEstimate)
        case 'grad'
            % Estimate noise from local gradients in the PSF image
            % This is assuming that the PSF is smooth
            [GxP,GyP]=gradient(P);
            % the median is required in order to get rid of the fast varying part
            % of the PSF (its core).
            Noise = median(sqrt( GxP(:).^2 + GyP(:).^2));
            Res(Isim).R0   = min(MatR(P<(InPar.NoiseFactor.*Noise)));
            if (~isempty(Res(Isim).R0))
                P(MatR>Res(Isim).R0) = 0;
            end
            Res(Isim).PSF = P./sum(P(:));
        case 'wmean'
            % Estimate PSF noise from weighted std
            NormCube = bsxfun(@times,Res(Isim).Cube(Res(Isim).Flag,:,:),1./Res(Isim).Flux(Res(Isim).Flag) );
            NormErr  = bsxfun(@plus,Res(Isim).Cube(Res(Isim).Flag,:,:),BackVariance(Res(Isim).Flag));
            NormErr  = bsxfun(@times,NormErr,1./Res(Isim).Flux(Res(Isim).Flag));
            [WM,WE]  = wmean(reshape(NormCube,Res(Isim).Ngood,StampWidth.^2),reshape(NormErr,Res(Isim).Ngood,StampWidth.^2),1);
            WM = reshape(WM,StampWidth,StampWidth);
            WE = reshape(WE,StampWidth,StampWidth)./sqrt(Res(Isim).Ngood);
            SN = WM./WE;
            Res(Isim).R0 = min(MatR(SN<(InPar.NoiseFactor)));
            if (~isempty(Res(Isim).R0))
                P(MatR>Res(Isim).R0) = 0;
            end
            Norm = sum(P(:));
            Res(Isim).PSF = P./Norm;
            Res(Isim).Ngood = sum(Res(Isim).Flag);
            Res(Isim).ErrPSF = WE./Norm;
        otherwise
            error('Unknown NoiseEstimate option');
    end
    
    % look for negative values
    Res.PSF(Res.PSF<0) = 0;
    
    if (InPar.Rad0<Inf)
        % Set PSF values outside some radius to 0
        Res(Isim).PSF(MatR>InPar.Rad0) = 0;
        Norm                           = sum(Res(Isim).PSF(:));
        Res(Isim).PSF                  = Res(Isim).PSF./Norm;
        
    end
    
    % count number of noz zero pixels in PSF
    Res(Isim).Nnonzero = sum(Res(Isim).PSF(:)>0);
    
    % Check number of non-zero pixels in PSF
    if (Res(Isim).Nnonzero<2)
        warning('PSF of images %d contains %d non-zero pixels\n',Isim,sum(Res(Isim).PSF>0)<2);
    end
    
    

    
    % Set output
    switch lower(InPar.OutType)
        case {'origsim','newsim','psf'}
            OutPSF(Isim).PSF    = Res(Isim).PSF;
            OutPSF(Isim).ErrPSF = Res(Isim).ErrPSF;
        otherwise
            error('Unknown OutType option');
    end
end



