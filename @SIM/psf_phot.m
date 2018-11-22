function [AstC,Col]=psf_phot(Sim,XY,varargin)
% PSF photometry for a list of sources position in a SIM object.
% Package: @SIM
% Description: PSF photometry for a list of sources position
%              on a single image (matrix).
%              The PSF photometry is done using a numerical PSF stamp.
%              For each source position, the PSF is sub pixel shifted to
%              the assumed source position and the flux and optionally
%              the background are fitted. The source position is not
%              fitted.
% Input  : - A SIM object. The image field should contain the background
%            component, and the BackIm and PSF fields should be available.
%            If the background field is empty then will attempt to fit the
%            background with the PSF fitting.
%          - PSF stamp must have odd size.
%            The PSF must be centered and on an whole pixel.
%            It is recomended that the PSF will be larger enough so
%            aperture corrections will be negligible. Note that the fit
%            itself may be done using only the core of the PSF (to avoid
%            influence by nearby sources).
%          - A list of positions in which to fit the PSF.
%            This is either a two column matrix of [X, Y],
%            or an AstCat object (which X/Y columns are defined by the
%            'ColNames' argument).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FitRad'   - Radius of the region in the PSF to use for the
%                         fit. Default is 2 pix.
%            'BackMethod'- Method by which to estimate background:
%                         'BackIm' - Take background from background image
%                                    in SIM. If doesn't exist then revert
%                                    to 'Ann'.
%                         'Ann'    - Estimate the background from an
%                                    annulus around the source.
%                         Default is 'Ann'.
%            'BackAnnRad' - Radii for background annulus [pix].
%                         Default is [7 9].
%            'BackAnnFun' - Function handle to use for background annulus
%                         estimation. Default is @median.
%            'BackSrc'  - A vector with background for each source.
%                         Default is empty.
%            'NormPSF'  - Normalize PSF to unity (true|false).
%                         This normalization is done on the whole PSF prior
%                         to selection of the PSF core (using 'FitRad') for
%                         fitting.
%                         Default is true.
%            'ColNames' - If the input sources position to measure is
%                         an AstCat object,
%                         then this is a cell array of column names
%                         that indicates which columns names contains
%                         the X/Y positions.
%                         Default is {'XWIN_IMAGE','YWIN_IMAGE'}.
%            'PreCalcPSF'- Pre calculate the shifted PSF in some finite
%                         grid (true). If false then will use the
%                         exact shift for each source.
%                         Default is false.
%            'Ngrid'    - The number of points in each axis in the grid
%                         in which to pre-calculate the shifted PSF.
%                         Default is 11.
%            'AlgoLSCOV'- Use lscov.m. Default is true.
%            'GetPSFPar'- Additional parameters to apass to
%                         ClassPSF/getpsf.m. Default is {}.
% Output : - A matrix or an AstCat object containing the resulted PSF
%            photometry catalog.
%          - Structure array containing the column names for the output
%            catalog.
% See also: psf_phot.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: SL=[135.1, 156.9 1000; 145.5 160.0 20000; 34.2 190.9 100000; 34.0 56.0 100000];
%          Im=image_art('StarList',SL,'Back',1000,'OutSim',false,'ParPSF',[1 1 0]); %,'AddNoiseSource',false,'AddNoiseBack',false);
%          [X,Y]=meshgrid((-5:1:5),(-5:1:5));
%          PSF = bivar_gauss(X,Y,[0 0 1 1 0 1 0]);
%          Res=psf_phot(Im,1000.*ones(size(Im)),PSF,SL(:,1:2));
%          SL = [rand(500,2).*200+10, ones(500,1).*1e4];
%          Im=image_art('StarList',SL,'Back',1000,'OutSim',false,'ParPSF',[1 1 0]); 
%          Res=psf_phot(Im,1000.*ones(size(Im)),PSF,SL(:,1:2));
%          [MatX,MatY]=meshgrid((20:20:200),(20:20:200));
%          SL = [MatX(:),MatY(:),1e4.*ones(size(MatX(:)))];
%          Im=image_art('StarList',SL,'Back',1000,'OutSim',false,'ParPSF',[3 3 0],'ImSize',[256 256]);
%          Res=psf_phot(Im,1000.*ones(size(Im)),PSF,SL(:,1:2));
% Reliable: 2
%--------------------------------------------------------------------------

ImageField = SIM.ImageField;
BackField  = SIM.BackField;

Nu  = 2;
T_P = tinv(normcdf(1,0,1),Nu);   % student t correction factor for 1/2 dof

CatField     = 'Cat';
ColField     = 'Col';

DefV.FitRad             = 2;
DefV.BackMethod         = 'Ann';
DefV.BackAnnRad         = [7 9];
DefV.BackAnnFun         = @median;
DefV.BackSrc            = [];
DefV.NormPSF            = true;
DefV.ColNames           = {'XWIN_IMAGE','YWIN_IMAGE'};
DefV.PreCalcPSF         = false;
DefV.Ngrid              = 11;
DefV.AlgoLSCOV          = false;
DefV.GetPSFPar          = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% AstCat output
AstC = AstCat(size(Sim));
            
Nsim = numel(Sim);
for Isim=1:1:Nsim
    Image = Sim(Isim).(ImageField);
    Back  = Sim(Isim).(BackField);
    PSF   = getpsf(Sim(Isim),InPar.GetPSFPar{:});
    
    
    
    ErrFactor = sqrt(sqrt(sum(PSF(:).^2)));

    % Get image size
    SizeIm = fliplr(size(Image));   % [X,Y] rather than [I,J]

    [MatX,MatY] = meshgrid((1:1:SizeIm(1)),(1:1:SizeIm(2)));
    
    % If XY is an AstCat object convert it into two column matrix [X,Y]
    if (AstCat.isastcat(XY))
        XY = col_get(XY,InPar.ColNames);
    end
    % shifts to nearest whole pixel
    RoundXY = round(XY);
    
    % make sure that PSF is normalized to 1
    if (InPar.NormPSF)
        PSF = PSF./sum(PSF(:));
    end

    % Perform the fit using only the PSF core
    SizeOrigPSF = size(PSF);
    Center = ceil(SizeOrigPSF.*0.5);
    PSF    = PSF(Center(1)-InPar.FitRad:1:Center(1)+InPar.FitRad,Center(2)-InPar.FitRad:1:Center(2)+InPar.FitRad);


    % half size of PSF - PSF is assumed to have odd-number size
    %SizeP    = size(PSF);
    SizePSF  = fliplr((size(PSF) - 1).*0.5);   % [X, Y] rather than [I,J]
    NumElPSF = numel(PSF);

    % select sources for which full PSF stamp is within image
    FlagSrc  = RoundXY(:,1)>SizePSF(1) & RoundXY(:,2)>SizePSF(2) & RoundXY(:,1)<(SizeIm(1)-SizePSF(1)) & RoundXY(:,2)<(SizeIm(2)-SizePSF(2));
    %RoundXY  = RoundXY(FlagSrc,:);
    %XY       = XY(FlagSrc,:);
    Delta    = XY - RoundXY;
    %if (~isempty(InPar.BackSrc))
    %    InPar.BackSrc = InPar.BackSrc(FlagSrc);
    %end

    
    if (isempty(Back))
        BackMethod = 'Ann';
    else
        BackMethod = InPar.BackMethod;
    end
    switch lower(BackMethod)
        case 'backim'
            % subtract background from image
            Image = Image - Back;
            % store the background for each source:
            BackXY = ImUtil.Im.image_posval(Back,RoundXY(:,1),RoundXY(:,2));
            
        otherwise
            % do nothing
    end
    

    % XY       = XY(FlagSrc,:);
    % RoundXY  = RoundXY(FlagSrc,:);
    % Delta    = Delta(FlagSrc,:);

    % number of sources for which to calculate PSF photometry:
    Nsrc = size(XY,1);

    % Initialize output catalog
    Col.X       = 1;
    Col.Y       = 2;
    Col.Flux    = 3;
    Col.FluxErr = 4;
    Col.Back    = 5;
    Col.BackErr = 6;
    Col.Chi2    = 7;
    Col.Chi2Back= 8;
    Col.Chi2CR  = 9;
    Cat = zeros(Nsrc,9).*NaN; % [X,Y,Flux,ErrFlux,Back,ErrBack,chi2]
    Cat(:,[Col.X, Col.Y]) = XY;


    % Construct catalog of shifted PSFs
    if (InPar.PreCalcPSF)
        Step = 1./(InPar.Ngrid-1);
        [ShX,ShY] = meshgrid((-0.5:Step:0.5),(-0.5:Step:0.5));
        SizeSh    = size(ShX);
        Nsh = numel(ShX);

        % for faster calculation - calc. NY,NX,Nr,Nc in first iteration
        Ish = 1;
        DX = ShX(Ish);
        DY = ShY(Ish);
        [CatPSF(Ish).PSF,NY,NX,Nr,Nc] = ImUtil.Im.image_shift_fft(PSF,DX,DY);
        for Ish=2:1:Nsh
            DX = ShX(Ish);
            DY = ShY(Ish);
            CatPSF(Ish).PSF = ImUtil.Im.image_shift_fft(PSF,DX,DY,NY,NX,Nr,Nc);
        end
    end


    % for each source in the input image
    FirstTime = true;
    for Isrc=1:1:Nsrc
        
        if (FlagSrc(Isrc))
            
            % for all sources withing bounderies
            % shift to nearest whole pixel
            DX = Delta(Isrc,1);
            DY = Delta(Isrc,2);

            % either use pre-constructed shifted PSF or shift PSF
            if (InPar.PreCalcPSF)
                % Use pre shifted PSF
                Ix = round((DX+0.5)./Step)+1;
                Iy = round((DY+0.5)./Step)+1;
                %Is = sub2ind(SizeSh,Iy,Ix);
                Is = Util.array.sub2ind_fast(SizeSh,Iy,Ix);
                ShiftedPSF = CatPSF(Is).PSF;
            else
                % shift PSF:
                % PSF does not contain bad pixels and therefore make sense
                % to shift PSF rather than the image
                if (FirstTime)
                    [ShiftedPSF,NY,NX,Nr,Nc] = ImUtil.Im.image_shift_fft(PSF,DX,DY);
                    FirstTime = false;
                else
                    [ShiftedPSF] = ImUtil.Im.image_shift_fft(PSF,DX,DY,NY,NX,Nr,Nc);
                end
            end

            % cut image
            % Note that sources near boundries are ignored
            %[CutImage] = ImUtil.Im.trim_image(Image,[RoundXY(Isrc,1),RoundXY(Isrc,2),SizePSF],'center',[]);
            [CutImage] = double(ImUtil.Im.trim_image(Image,[RoundXY(Isrc,1),RoundXY(Isrc,2),SizePSF],'center',[])); % Na'ama, 20180828

            % Index of center of CutImage
            IndCenterIm = ceil(numel(CutImage).*0.5);

            %---------------
            %--- Fit PSF ---
            %---------------
            switch lower(BackMethod)
                case 'backim'
                    % use the background image
                    % background is already subtracted
                    % fit only the flux
                    BackEst = BackXY(Isrc);
                    BackErrEst = NaN;
                case 'ann'
                    % estimate background from annulus around source
                    if (isempty(InPar.BackSrc))
                        % estimate the background for each source
                        % old - slow
%                         MatR2 = (MatX - XY(Isrc,1)).^2 + (MatY - XY(Isrc,2)).^2;
%                         FlagBack = MatR2>InPar.BackAnnRad(1).^2 & MatR2<=InPar.BackAnnRad(2).^2;
%                         BackEst  = InPar.BackAnnFun(Image(FlagBack));
%                         BackErrEst = std(Image(FlagBack))./sqrt(sum(FlagBack(:)));
                         
                        % fast estimate
                        [BackEst,BackErrEst] = ImUtil.Im.annulus_properties(Image,XY(Isrc,1),XY(Isrc,2),InPar.BackAnnRad,InPar.BackAnnFun,@std,MatX,MatY);
                    else
                        BackEst    = InPar.BackSrc(Isrc);
                        BackErrEst = NaN;
                    end

                    CutImage = CutImage - BackEst;
                otherwise
                    error('Unknown BackMethod option');
            end

            % ShiftedPSF(:) is the design matrix
            Par    = ShiftedPSF(:)\CutImage(:);
            Resid  = CutImage(:) - ShiftedPSF(:)*Par;
            ParErr = std(Resid);  %./sqrt(numel(ShiftedPSF));


            Cat(Isrc,Col.Flux)    = Par;
            Cat(Isrc,Col.FluxErr) = ParErr;
            Cat(Isrc,Col.Back)    = BackEst;
            Cat(Isrc,Col.BackErr) = BackErrEst;


            % \chi^2 for point source (PSF) model
            % note that CutImage is background subtracted
            Cat(Isrc,Col.Chi2) = sum(Resid.^2./(CutImage(:)+BackEst));
            % \chi^2 for no source (i.e., background)
            Cat(Isrc,Col.Chi2Back) = sum((CutImage(:)-BackEst).^2./(CutImage(:)+BackEst));
            % \chi^2 for delta function (CR) model
            Cat(Isrc,Col.Chi2CR) = Cat(Isrc,Col.Chi2Back) - CutImage(IndCenterIm).^2./(CutImage(IndCenterIm)+BackEst);

        end
    end

    AstC(Isim).(CatField) = Cat;
    AstC(Isim).(ColField) = Col;

end
AstC = col2colcell(AstC);