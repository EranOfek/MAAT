function [CPSF,Res]=psf_combiner(CubePSF,varargin)
% Combine a cube of PSFs
% Package: ImUtil.Im
% Description: Given a cube of PSF stamps, align the PSFs, remove bad PSFs
%              and combine them into a single PSF.
% Input  : - A Cube of PSFs where the first index is the image index.
%            The PSF needs to be roughly in the center of each stamp.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'BackVec' - A vecor of background value per each PSF stamp.
%                        If 0, then equivalent to background already
%                        subtracted. If empty, then estimate background
%                        from outer region of stamp.
%            'BackAnnrad' - The outer radius of the annulus to use for
%                        background estimation in each stamp.
%                        If empty, then set outer radius to stamp radius.
%                        Default is empty.
%            'BackAnnWidth' - The width of the background annulus [pix].
%                        Default is 2.
%            'MomRadius' - Radius in which to calculate moments.
%                        If empty then set to BackAnnRad-BackAnnWidth.
%                        Default is empty.
%            'MomSigma' - Moment calculation weight sigma [pix].
%                        Default is 1.5.
%            'MomMaxIter' - Maximum number of iterations for moments
%                        calculation. Default is 5.
%            'AdjustMomSigma' - If true, then estimate Sigma from PSF width
%                        and rerun the moment calculation with the
%                        estimated sigma. Default is true.
%            'BackFun'   - Function handle to use for the background
%                        calculation in the annulus.
%                        Default is @median.
%            'AlignPSFs' - Align PSF using fft sub-pixel shifts.
%                        Default is true.
%            'CleanMoments' - A logical flag indicating if to remove PSF
%                        stamps which 2nd moments are outliers.
%                        Default is true.
%            'CleanMomentsNsig' - Number of sigma, for clean moments sigma
%                        clipping. Default is 3.
%            'CleanMultiPeak' - Remove stars with stamps include bright
%                        secondary local maxima, which S/N ratio is larger
%                        than 'MultiPeakThresh'.
%            'ConnMax' - Local maxima connectivity parameter for
%                        imregionalmax. Default is 8.
%            'MultiPeakThresh' -  If more than one source with maximum
%                        S/N per pixel larger than this value exist in
%                        the source PSF stamp than remove stamp.
%                        Default is 3.
%            'CleanBack' - A logical flag indicating if to remove PSF
%                        stamps which have bad background level or std.
%                        Default is true.
%            'Backquantile' - Remove stamps with background level above
%                        this quantile level. Default is 0.9.
%            'BackStdquantile' - emove stamps with background std above
%                        this quantile level. Default is 0.9.    
%            'CleanResid' - A logical flag indicating if to remove PSFs
%                        with bad residuals based on total chi^2 of each
%                        PSF relative to the mean PSF and maximum deviation
%                        pixel in each PSF.
%                        Default is true.
%            'Chi2quantile' - Upper quantile of \chi^2 values above to
%                        remove bad PSFs. Default is 0.9.
%            'MaxDquantile' - Upper quantile of maximum deviating pixels
%                        in std units above which to remove the PSF stamp.
%                        Default is 0.9.
%            'Method'   - Function handle for PSF combining method.
%                        Default is wmean.
%            'FinalBackSub' - A logical flag indicating if to subtract
%                        background residuals from final PSF.
%                        Default is true.
%            'MinSN'   - S/N of the final PSF below which to set the PSF
%                        value to zero. If empty, then do nothing.
%                        Default is 3.
%            'ErrMethod' - For Method='wmean' this set the error estimation
%                        method.
%                        'err' | 'stdn'. Default is 'stdn'.
%            'NsigmaRad0' - Any pixel in the PSF which distance from the
%                        PSF center is larger than the second moment sigma
%                        multiplied by this factor is set to 0.
%                        If empty, do nothing. Default is 5.
%            'SetNegative' - A logixcal flag indicating if to set negative
%                        PSF values to 0. Default is true.
%            'Norm'    - Normalize PSF, such that its sum will equal this
%                        value. If empty, do nothing. Default is 1.
%            'Verbose' - Verbose. Default is false.
% Output : - A class PSF object with the combined PSF
%          - A structure array of summary information.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: F=rand(100,1).*50+50;
%          for I=1:1:100, Cube(I,:,:)=poissrnd(Kernel2.gauss.*F(I)); end
%          P=ImUtil.Im.psf_combiner(Cube)
% Reliable: 2
%--------------------------------------------------------------------------


DefV.BackVec              = [];  % if empty, then estimate background - 0 for background already subtracted.
DefV.BackAnnRad           = [];  % if empty choose half stamp size (outer radius)
DefV.BackAnnWidth         = 2;
DefV.MomRadius            = [];  % default is BackAnnRad-BackAnnWidth
DefV.MomSigma             = 1.5; % [pix]
DefV.MomMaxIter           = 5;
DefV.AdjustMomSigma       = true;
DefV.BackFun              = @median;
DefV.AlignPSFs            = true;
DefV.CleanMoments         = true;
DefV.CleanMomentsNsig     = 3;
DefV.CleanMultiPeak       = false;
DefV.ConnMax              = 8;
DefV.MultiPeakThresh      = 3;
DefV.CleanBack            = true;
DefV.Backquantile         = 0.9;
DefV.BackStdquantile      = 0.9;
DefV.CleanResid           = true;
DefV.Chi2quantile         = 0.9;
DefV.MaxDquantile         = 0.9;
DefV.Method               = 'wmean'; %@mean; %'wmedian';  %@mean;
DefV.FinalBackSub         = true;
DefV.MinSN                = 3;
DefV.ErrMethod            = 'stdn';
DefV.SetNegative          = true;
DefV.NsigmaRad0           = 6; % 5; %3;
DefV.Norm                 = 1;
DefV.Verbose              = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


SizePSF = size(CubePSF);
Npsf    = SizePSF(1);
SizePSF = SizePSF(2:3);

Center = (SizePSF+1).*0.5;
[MatX,MatY] = meshgrid((1:1:SizePSF(2))-Center(2),(1:1:SizePSF(1))-Center(1));

if (isempty(InPar.BackAnnRad))
    InPar.BackAnnRad = min(floor(SizePSF.*0.5));
end

if (isempty(InPar.MomRadius))
    InPar.MomRadius = InPar.BackAnnRad - InPar.BackAnnWidth;
end



% Estimate the background for each PSF
if isempty(InPar.BackVec)
    % Estimate the background
    MatR2 = MatX.^2 + MatY.^2;
    
    % select pixels in the background annulus
    FlagBack = MatR2>(InPar.BackAnnRad-InPar.BackAnnWidth).^2 & MatR2<=InPar.BackAnnRad.^2;
    FlagBack = FlagBack(:);
    InPar.BackVec = zeros(Npsf,1);
    BackStd       = zeros(Npsf,1);
    for Ipsf=1:1:Npsf
        % calc background and std of background in each PSF candidate
        CubePSF1 = squeeze(CubePSF(Ipsf,:,:));
        InPar.BackVec(Ipsf) = InPar.BackFun(CubePSF1(FlagBack));
        BackStd(Ipsf)       = std(CubePSF1(FlagBack));
    end
end

% in case that background is scalar: duplicate
if (numel(InPar.BackVec)==1)
    InPar.BackVec = InPar.BackVec.*ones(Npsf,1);
end

% Subtract the background for each PSF
CubePSFbs = zeros(Npsf,SizePSF(1),SizePSF(2));  % back. sub. PSF cube
for Ipsf=1:1:Npsf
    CubePSFbs(Ipsf,:,:) = CubePSF(Ipsf,:,:) - InPar.BackVec(Ipsf);
end


% Calculate 1st and 2nd moments
for Ipsf=1:1:Npsf
    [Mom(Ipsf),Mom2(Ipsf)] = ImUtil.Im.im_moments(squeeze(CubePSFbs(Ipsf,:,:)),Center(2),Center(1),InPar.MomRadius,InPar.MomSigma,InPar.MomMaxIter);
end

% second iteration for 1st and 2nd moment using estimated sigma
% The reason for this is that im_moments uses a Gaussian weight function.
% This step allows you to adjust the weight function width based on the
% width found in the first iteration
if (InPar.AdjustMomSigma)
    InPar.MomSigma = sqrt(nanmedian([Mom2.X2]).^2 + nanmedian([Mom2.Y2]).^2);
    % Calculate 1st and 2nd moments
    for Ipsf=1:1:Npsf
        [Mom(Ipsf),Mom2(Ipsf)] = ImUtil.Im.im_moments(squeeze(CubePSFbs(Ipsf,:,:)),Center(2),Center(1),InPar.MomRadius,InPar.MomSigma,InPar.MomMaxIter);
    end
end

% align the PSFs
if (InPar.AlignPSFs)
    for Ipsf=1:1:Npsf
        DX = Mom(Ipsf).X - Center(2);
        DY = Mom(Ipsf).Y - Center(1);
        %CubePSFbs(Ipsf,:,:) = ImUtil.Im.image_shift_fft(squeeze(CubePSFbs(Ipsf,:,:)),DX,DY);
        if (Ipsf==1)
           [CubePSFbs(Ipsf,:,:),NY,NX,Nr,Nc] = ImUtil.Im.image_shift_fft(squeeze(CubePSFbs(Ipsf,:,:)),DX,DY);
        else
           [CubePSFbs(Ipsf,:,:)] = ImUtil.Im.image_shift_fft(squeeze(CubePSFbs(Ipsf,:,:)),DX,DY,NY,NX,Nr,Nc);
        end
        CubePSFbs(Ipsf,:,:) = CubePSFbs(Ipsf,:,:);
    end
end

%----------------------
%--- clean bad PSFs ---
%----------------------

% clean by moments
if (InPar.CleanMoments)
    MomX1 = [Mom.X].';
    MomY1 = [Mom.Y].';
    MomX2 = [Mom2.X2].';
    MomY2 = [Mom2.Y2].';
    MomXY = [Mom2.XY].';
    
    StdX1 = Util.stat.rstd(MomX1);
    StdY1 = Util.stat.rstd(MomY1);
    
    MedX2 = nanmedian(MomX2);
    StdX2 = Util.stat.rstd(MomX2);
    MedY2 = nanmedian(MomY2);
    StdY2 = Util.stat.rstd(MomY2);
    MedXY = nanmedian(MomXY);
    StdXY = Util.stat.rstd(MomXY);
    
    % Flag good stars
    FlagMom1 = abs(Center(2)-MomX1)<(StdX1.*InPar.CleanMomentsNsig) & ...
               abs(Center(1)-MomY1)<(StdY1.*InPar.CleanMomentsNsig);
    FlagMom2 = abs(MomX2-MedX2)<(StdX2.*InPar.CleanMomentsNsig) & ...
               abs(MomY2-MedY2)<(StdY2.*InPar.CleanMomentsNsig) & ...
               abs(MomXY-MedXY)<(StdXY.*InPar.CleanMomentsNsig); 
    
    FlagMom = FlagMom1 & FlagMom2;
    Npsf = sum(FlagMom);
    if (InPar.Verbose)
        fprintf('PSF stars with good moments: %d\n',Npsf);
    end
    
    CubePSFbs     = CubePSFbs(FlagMom,:,:);
    InPar.BackVec = InPar.BackVec(FlagMom);
    BackStd       = BackStd(FlagMom);
    
    Mom           = Mom(FlagMom);
    Mom2          = Mom2(FlagMom);
end

% clean by multiple peaks flux ratio
FlagMultiPeak = true(Npsf,1);
if (InPar.CleanMultiPeak)
    for Ipsf=1:1:Npsf
        Tmp = squeeze(CubePSFbs(Ipsf,:,:));
        MaxBW = imregionalmax(Tmp,InPar.ConnMax);
        if (sum(sort(Tmp(MaxBW)./std(Tmp(:)))>InPar.MultiPeakThresh)>1)
            % more than one peak found - bad source
            FlagMultiPeak(Ipsf) = false;
        
        end
        
    end
    CubePSFbs     = CubePSFbs(FlagMultiPeak,:,:);
    Npsf          = sum(FlagMultiPeak);
    Mom           = Mom(FlagMultiPeak);
    Mom2          = Mom2(FlagMultiPeak);
    InPar.BackVec = InPar.BackVec(FlagMultiPeak);
    BackStd       = BackStd(FlagMultiPeak);
end



% clean by background level and std
if (InPar.CleanBack)
    Res.MaxBackLevel = quantile(InPar.BackVec,InPar.Backquantile);
    FlagBackLevel    = InPar.BackVec<=Res.MaxBackLevel;
    
    Res.MaxBackStd   = quantile(BackStd,InPar.BackStdquantile);
    FlagBackStd      = BackStd<=Res.MaxBackStd;
    
    FlagBack         = FlagBackLevel & FlagBackStd;
    
    Npsf          = sum(FlagBack);
    CubePSFbs     = CubePSFbs(FlagBack,:,:);
    InPar.BackVec = InPar.BackVec(FlagBack);
    BackStd       = BackStd(FlagBack);
    Mom           = Mom(FlagBack);
    Mom2          = Mom2(FlagBack);
    
end


% clean by residuals
if (InPar.CleanResid)
    % rough flux estimate
    Flux = zeros(Npsf,1);
    for Ipsf=1:1:Npsf
        Flux(Ipsf) = sum(sum(CubePSFbs(Ipsf,:,:)));
    end
    % renormalize PSF cube such the sum of each PSF is 1
    NormCubePSFbs = bsxfun(@times,CubePSFbs,1./Flux);
    
    % get rough PSF using median
    PSF     = nanmedian(NormCubePSFbs,1);
    % get robust estimation of noise in PSF
    StdPSF  = Util.stat.rstd(NormCubePSFbs,1);
    
    % residual for each PSF
    Resid = bsxfun(@minus,NormCubePSFbs,PSF);
    % The residuals divided by the noise
    ResidStd = bsxfun(@times,Resid,1./StdPSF);
    % \chi^2
    Chi2       = sum(sum(ResidStd.^2,3),2);
    Res.Chi2qunatile = quantile(Chi2,InPar.Chi2quantile);
    FlagChi2   = Chi2<Res.Chi2qunatile;

    % search for individual stars with large deviations
    MaxDev = zeros(Npsf,1);
    for Ipsf=1:1:Npsf
        MaxDev(Ipsf) = max(max(abs(ResidStd(Ipsf,:,:))));
    end
    Res.MaxDqunatile = quantile(MaxDev,InPar.MaxDquantile);
    FlagMaxD   = MaxDev<=Res.MaxDqunatile;
    
    
    FlagResid     = FlagChi2 & FlagMaxD;
    Npsf          = sum(FlagResid);
    CubePSFbs     = CubePSFbs(FlagResid,:,:);
    InPar.BackVec = InPar.BackVec(FlagResid);
    BackStd       = BackStd(FlagResid);
    Mom           = Mom(FlagResid);
    Mom2          = Mom2(FlagResid);
    
    if (InPar.Verbose)
        fprintf('PSF stars with good residuals: %d\n',Npsf);
    end
    
end

% Estimated flux per PSF
Flux = zeros(Npsf,1);
for Ipsf=1:1:Npsf
    Flux(Ipsf) = sum(sum(CubePSFbs(Ipsf,:,:)));
end
NormCubePSFbs = bsxfun(@times,CubePSFbs,1./Flux);

%------------------------
%--- Combine the PSFs ---
%------------------------
if (isa(InPar.Method,'function_handle'))
    % InPar.Method is a function that get a cube and sum over the first
    % dimension
    PSF = squeeze(InPar.Method(NormCubePSFbs));
    %ErrPSF = squeeze(std(NormCubePSFbs))./sqrt(Npsf);
    ErrPSF = squeeze(Util.stat.rstd(NormCubePSFbs))./sqrt(Npsf);
else
    switch lower(InPar.Method)
        case 'wmedian'
            % Construct PSF using weighted median
            PSF = zeros(SizePSF);
            for Ix=1:1:SizePSF(2)
                for Iy=1:1:SizePSF(1)
                    Val = squeeze(NormCubePSFbs(:,Iy,Ix));
                    Err = sqrt(abs(squeeze(CubePSFbs(:,Iy,Ix)) + InPar.BackVec));
                    PSF(Iy,Ix) = Util.stat.wmedian(Val, Err./Flux);
                   
                end
            end
            ErrPSF = nan(size(PSF));
%         case 'wmean'
%             % Construct PSF using weighted mean
%             PSF = zeros(SizePSF);
%             for Ix=1:1:SizePSF(2)
%                 for Iy=1:1:SizePSF(1)
%                     Val = squeeze(CubePSFbs(:,Iy,Ix));
%                     Err = sqrt(abs(squeeze(NormCubePSFbs(:,Iy,Ix)) + InPar.BackVec));
%                     PSF(Iy,Ix) = Util.stat.wmean(Val, Err./Flux);
%                 end
%             end
%             ErrPSF = nan(size(PSF));
        case 'lsq'
            % Least Square fitting
            % Under development
            
            CubePSFbs
            % Pij - The PSF
            % i - image
            % j - pixel
            Npix = prod(SizePSF);
            Pij = reshape(CubePSFbs,[Npsf Npix]);
            
            %Pij=[1 2 3;4 5 6;7 8 9;10 11 12]
            %[Npsf,Npix] = size(Pij);
            
            Hdiag  = zeros(Npsf.*Npix,Npsf);
            Hblock = zeros(Npsf.*Npix,Npix); 
            
            Hdiag  = sparse([],[],[],Npsf.*Npix,Npsf,Npsf.*Npix);
            Hblock = sparse([],[],[],Npsf.*Npix,Npix,Npsf.*Npix);
            for Iblock=1:1:Npix
                I1 = (Iblock-1).*Npsf+1;
                I2 = Iblock.*Npsf;
                Hdiag(I1:I2,:)       = diag(Pij(:,Iblock));
                Hblock(I1:I2,Iblock) = 1;  
            end
            H = [Hdiag, Hblock];
            eigs(H*H',1,'SM')
            
        case 'wmean'
            
            % no squeeze:
            if (size(NormCubePSFbs,1)==1)
                % Cube contains a single PSF
                PSF    = squeeze(NormCubePSFbs);
                ErrPSF = NaN.*PSF;
                warning('PSF is based in a single star');
            else
                ApproxPSF       = squeeze(median(NormCubePSFbs));  % for weight
                ApproxPSF       = ApproxPSF./sum(ApproxPSF(:));
                %NormFlux        = zeros(Npsf,1);

                Npix = prod(SizePSF);
                % i - pixel
                % j - image
                Pij = reshape(CubePSFbs,[Npsf Npix]).';
                FlagFit = MatR2<2.5.^2;
                % PSF fitting of inner core
                ApproxPSF1     = ApproxPSF;
                Niter = 4;
                for Iiter=1:1:Niter
                    NormFlux = (ApproxPSF1(FlagFit)\Pij(FlagFit,:)).';
                    % normalize to ~unity
                    NormCubePSFbs1 = bsxfun(@times,CubePSFbs,1./NormFlux);
                    % median combine normalized PSFs
                    %ApproxPSF1       = median(NormCubePSFbs1);  % for weight
                    Err = 1./sqrt(repmat(NormFlux+InPar.BackVec,1,SizePSF(1),SizePSF(2)));
                    [ApproxPSF1,ErrPSF,StdPSF] = Util.stat.wmean(NormCubePSFbs1,Err,1,false);

                    % renormalize to 1
                    ApproxPSF1       = ApproxPSF1./sum(ApproxPSF1(:));
                end
                PSF = squeeze(ApproxPSF1);
                switch lower(InPar.ErrMethod)
                    case 'err'
                        ErrPSF = squeeze(ErrPSF);
                    case 'stdn'
                        % It is important to use this option when the
                        % original images are background subtracted.
                        % In this case the noise properties can not bes
                        % estimated from the background value...
                        ErrPSF = squeeze(StdPSF)./sqrt(Npsf);
                    otherwise
                        error('Unknown ErrMethod option');
                end
            end
        otherwise
            error('Unknown Method option');
    end
end

% subtract and background residuals
MatR2 = MatX.^2 + MatY.^2;
    
FlagBack       = MatR2>(InPar.BackAnnRad-InPar.BackAnnWidth).^2 & MatR2<=InPar.BackAnnRad.^2;
FlagBack       = FlagBack(:);
Res.BackPSF    = InPar.BackFun(PSF(FlagBack));
Res.StdBackPSF = std(CubePSF1(FlagBack));
if (InPar.FinalBackSub)
    PSF = PSF - Res.BackPSF;
end

% set low S/N values to zero
if (~isempty(InPar.MinSN))
    SN = PSF./ErrPSF;
    PSF(~isnan(SN) & SN<InPar.MinSN) = 0;
end

% set negative values to zero
if (InPar.SetNegative)
    PSF(PSF<0) = 0;
end

[Res.Mom,Res.Mom2] = ImUtil.Im.im_moments(PSF,Center(2),Center(1),InPar.MomRadius,InPar.MomSigma,InPar.MomMaxIter);
Res.SigmaPSF = sqrt(median([Mom2.X2]).^2 + median([Mom2.Y2]).^2);

Res.Npsf = Npsf;

if (~isempty(InPar.NsigmaRad0))
    PSF(MatR2>(InPar.NsigmaRad0.*Res.SigmaPSF).^2) = 0;
end
    
% Normalize PSF
if (~isempty(InPar.Norm))
    Norm   = sum(PSF(:));
    PSF    = PSF.*InPar.Norm./Norm;
    ErrPSF = ErrPSF.*InPar.Norm./Norm;
end

% prepare output
CPSF        = ClassPSF;
CPSF.PSF    = PSF;
CPSF.ErrPSF = ErrPSF;

if (Npsf==0)
    error('Failed to combine PSF');
end