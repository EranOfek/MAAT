function [Res]=fit_transform(Ref,Cat,TranC,varargin)
% Fit astrometric transformation
% Package: +ImUtil/+pattern
% Description: 
% Input  : - Reference catalog in a matrix format.
%            By default the X/Y (RA/Dec) coordinates are in columns 1,2,
%            and the magnitude and color are in column 5,6 (See
%            ColRefX,...).
%          - Catalog to fit to the reference in a matrix format.
%          - A transformation class (TranClass) object that defines the
%            fitted transformation.
%            Default is 'affine_tt_cheby2_4'.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColRefX' - Column of X coordinate in Ref catalog.
%                        Default is 1.
%            'ColRefY' - Column of Y coordinate in Ref catalog.
%                        Default is 2.
%            'ColRefMag' - Column of magnitude in Ref catalog.
%                        Default is 5.
%            'ColRefColor' - Column of color in Ref catalog.
%                        Default is 6.
%            'ColCatX' - Column of X coordinate in catalog.
%                        Default is 1.
%            'ColCatY' - Column of Y coordinate in catalog.
%                        Default is 2.
%            'PixScale' - Pixel scale in arcsec/pix. Used for the rms
%                        estimate in the image sub blocks. 
%            'CooUnits' - Coordinate units. Default is 'deg'.
%            'RefSelectFun' - Function handle to apply on the reference
%                        catalog for selection of sources that will
%                        participate in the fit.
%                        Default is @(R,Par) R(:,Par)<19 & R(:,Par)>12;
%            'RefSelectFunPar' - Cell array of parameters to pass to the
%                        'RefSelectFun'. Default is {5}.
%            'NormXY' - Normalization of the X/Y coordinates prior to the
%                       fit. Default si 4096.
%            'PolyMagDeg' - Max polynomial degree for magnitude fit.
%                       Default is 3.
%            'StepMag' - Magnitude bin size in rms vs. mag diagaram.
%                      Default is 0.1 mag.
%                      This is being used in order to estimate the rms as a
%                      function of magnitude, which in turn is used as the
%                      error estimate on the measurment in the second
%                      iteration.
%            'Niter' - Number of extra iterations to perform after the
%                      first fit. Default is 1.
%            'SigClip' - Sigma clipping to apply after the basic iteration.
%                      Default is 5.
%            'MaxResid' - Maximum allowed residuals in pixels.
%                      After the basic fit, larger residuals will be
%                      removed.
%                      Default is 0.1 pix.
%            'Plot' -  Default is false.
%            'ImSize' -  Default is [].
%            'BlockSize' - Default is [].
%            'BufferSize' Default is 10.
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=ImUtil.pattern.fit_transform(Ref,Cat)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3)
    TranC = [];
end
if (isempty(TranC))
    TranC = 'affine_tt_cheby2_4';
end


DefV.ColRefX              = 1;
DefV.ColRefY              = 2;
DefV.ColRefMag            = 5;
DefV.ColRefColor          = 6;
DefV.ColCatX              = 1;
DefV.ColCatY              = 2;
DefV.PixScale             = 1;   % arcsec/pix
DefV.CooUnits             = 'deg';
DefV.RefSelectFun         = @(R,Par) R(:,Par)<19 & R(:,Par)>12;
DefV.RefSelectFunPar      = {5};
DefV.NormXY               = 4096;
DefV.PolyMagDeg           = 3;
DefV.StepMag              = 0.1;
DefV.Niter                = 1;
DefV.SigClip              = 5;
DefV.MaxResid             = 1./3600;  % assuming coordinates are in deg
DefV.Plot                 = false;
DefV.ImSize               = []; %[2048 4096];  %X/Y
DefV.BlockSize            = []; %[512 512];    %X/Y
DefV.BufferSize           = 10;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% select sources in RefCat
if (~isempty(InPar.RefSelectFun))
    Flag = InPar.RefSelectFun(Ref,InPar.RefSelectFunPar{:});
    Ref  = Ref(Flag,:);
    Cat  = Cat(Flag,:);
    
else
    Flag = true(size(Ref,1),1);    
end

RefX = Ref(:,InPar.ColRefX);
RefY = Ref(:,InPar.ColRefY);
CatX = Cat(:,InPar.ColCatX);
CatY = Cat(:,InPar.ColCatY);
if (isempty(InPar.ColRefMag))
    RefMag = [];
else
    RefMag = Ref(:,InPar.ColRefMag);
end
if (isempty(InPar.ColRefColor))
    RefColor = [];
else
    RefColor = Ref(:,InPar.ColRefColor);
end

N = numel(RefX);


% Deal with TranC
if (isa(TranC,'TranClass'))
    
        
elseif (isnumeric(TranC))
    
elseif (ischar(TranC))
    switch lower(TranC)
        case 'shift'
            TranC = TranClass({@FunOne},...
                              {@FunOne});

        case 'affine'
            % need to use FunM*
            TranC = TranClass({@FunOne, [],@FunX,[],@FunY,[]},...
                              {@FunOne, [],@FunX,[],@FunY,[]});

        case 'affine_tt'
            TranC = TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunTiltXp,[],@FunTiltXn,[]},...
                              {@FunOne, [],@FunX,[],@FunY,[],@FunTiltYp,[],@FunTiltYn,[]});

        case 'affine_tt_cheby2_3'
            ChebyOrder = [2 2;  3 0; 3 1; 3 2; 0 3; 1 3; 2 3; 3 3]';
            TranC = TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunTiltXp,[],@FunTiltXn,[],@FunPolyChebyshev2XY,ChebyOrder},...
                              {@FunOne, [],@FunX,[],@FunY,[],@FunTiltYp,[],@FunTiltYn,[],@FunPolyChebyshev2XY,ChebyOrder });
                  
        case 'affine_tt_cheby2_4'
            ChebyOrder = [2 2;  3 0; 3 1; 3 2; 0 3; 1 3; 2 3; 3 3; 4 0; 0 4; 4 1; 1 4; 4 2; 2 4; 4 3; 3 4; 4 4]';
            TranC = TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunTiltXp,[],@FunTiltXn,[],@FunPolyChebyshev2XY,ChebyOrder},...
                              {@FunOne, [],@FunX,[],@FunY,[],@FunTiltYp,[],@FunTiltYn,[],@FunPolyChebyshev2XY,ChebyOrder });
                          
        case 'affine_tt_cheby2_5'
            ChebyOrder = [2 2;  3 0; 3 1; 3 2; 0 3; 1 3; 2 3; 3 3; 4 0; 0 4; 4 1; 1 4; 4 2; 2 4; 4 3; 3 4; 4 4; 1 5; 5 1; 2 5; 5 2; 3 5; 5 3; 4 5; 5 4; 5 5]';
            TranC = TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunTiltXp,[],@FunTiltXn,[],@FunPolyChebyshev2XY,ChebyOrder},...
                              {@FunOne, [],@FunX,[],@FunY,[],@FunTiltYp,[],@FunTiltYn,[],@FunPolyChebyshev2XY,ChebyOrder });                  

        otherwise
            error('Unknown TranC option');
    end
else
    error('Unknown TranC option');
end


Res.TranC = TranC;

%H = design_matrix(TranC,'X',RefX./InPar.NormXY,'Y',RefY./InPar.NormXY);
%---
H = design_matrix(TranC,'X',CatX./InPar.NormXY,'Y',CatY./InPar.NormXY);

Ix = 1;
Iy = 2;

%ParX = H{Ix}\CatX;
%ParY = H{Iy}\CatY;
%---
ParX = H{Ix}\RefX;
ParY = H{Iy}\RefY;

%ResidX = CatX - H{1}*ParX;
%ResidY = CatY - H{2}*ParY;
%---
ResidX = RefX - H{1}*ParX;
ResidY = RefY - H{2}*ParY;

Resid  = sqrt(ResidX.^2 + ResidY.^2);

Res.TranC  = TranC;

Res.rrms1x = Util.stat.rstd(ResidX);
Res.rrms1y = Util.stat.rstd(ResidY);
Res.rrms1  = Util.stat.rstd(Resid);
Res.Nsrc1  = numel(Resid);

if (~isempty(RefMag))
    % Magnitude is available
    % fit residuals as a function of magnitude
    [Par,SPar] = polyfit(RefMag,log10(Resid),InPar.PolyMagDeg);
    
    FlagG = true(numel(RefMag),1);

    for Iiter=1:1:InPar.Niter
        [Par,SPar] = polyfit(RefMag(FlagG),log10(Resid(FlagG)),InPar.PolyMagDeg);
        
        ErrY     = SPar.normr^2/SPar.df;
        ErrYmag  = 10.^polyval(Par,RefMag);
        Xmag  = (min(RefMag):InPar.StepMag:max(RefMag))';
        ErrYmagX = 10.^polyval(Par,Xmag);
        [MinErr,MinI] = min(ErrYmagX);
        %Xmag(MinI)
        FlagG =  Resid<(ErrYmag + ErrYmag.*ErrY.*InPar.SigClip) & Resid<(InPar.MaxResid);
        
        %[ParX, ParErrX] = lscov(H{Ix}(FlagG,:), CatX(FlagG),  1./(ErrYmag(FlagG).^2));
        %[ParY, ParErrY] = lscov(H{Iy}(FlagG,:), CatY(FlagG),  1./(ErrYmag(FlagG).^2));
        %---
        [ParX, ParErrX] = lscov(H{Ix}(FlagG,:), RefX(FlagG),  1./(ErrYmag(FlagG).^2));
        [ParY, ParErrY] = lscov(H{Iy}(FlagG,:), RefY(FlagG),  1./(ErrYmag(FlagG).^2));
        

        %ResidX = CatX - H{1}*ParX;
        %ResidY = CatY - H{2}*ParY;
        %---
        ResidX = RefX - H{1}*ParX;
        ResidY = RefY - H{2}*ParY;
        
        Resid  = sqrt(ResidX.^2 + ResidY.^2);

       
    end
    [Par,SPar] = polyfit(RefMag(FlagG),log10(Resid(FlagG)),InPar.PolyMagDeg);
        
    %ResFitNoise = Util.fit.fit_noise_model(Resid(FlagG),RefMag(FlagG))
    
    ErrY     = SPar.normr^2/SPar.df;
    ErrYmag  = 10.^polyval(Par,RefMag);
    Xmag  = (min(RefMag(FlagG)):InPar.StepMag:max(RefMag(FlagG)))';
    ErrYmagX = 10.^polyval(Par,Xmag);
    [MinErr,MinI] = min(ErrYmagX);
    
    if (any(isnan(Par)) || any (isinf(Par)))
        Imin = [];
    else
        RootsP = abs(roots([3.*Par(1), 2.*Par(2), Par(3)]));
        Par2   = [6.*Par(1), 2.*Par(2)];
        Imin   = find(polyval(Par2,RootsP)>0);
    end

    if (~isempty(Imin))
        MagMin = RootsP(Imin(1));

        if (isempty(MagMin))
            LocalMinErr = NaN;
        else
            if (MagMin>min(RefMag) && MagMin<max(RefMag))
                LocalMinErr = 10.^polyval(Par,MagMin);
            else
                LocalMinErr = NaN;
            end
        end
    else
        LocalMinErr = NaN;
    end
        
    if (InPar.Plot)
        figure;
        semilogy(RefMag(FlagG),(Resid(FlagG)),'.')    
        hold on
        plot(Xmag,10.^polyval(Par,Xmag),'k-') 
    end

    Res.ResidX = ResidX;
    Res.ResidY = ResidY;
    Res.Resid  = Resid;
    Res.RefMag = RefMag;
    Res.RefColor = RefColor;
    Res.CatX   = CatX;
    Res.CatY   = CatY;
    Res.RefX   = RefX;
    Res.RefY   = RefY;
    Res.FlagG  = FlagG;
    Res.rmsN   = std(Resid(FlagG));
    Res.rrmsNx = Util.stat.rstd(ResidX(FlagG));
    Res.rrmsNy = Util.stat.rstd(ResidY(FlagG));
    Res.rrmsN  = Util.stat.rstd(Resid(FlagG));
    Res.NsrcN  = sum(FlagG);
    Res.wmedrms  = NaN; %Util.stat.wmedian(Resid(FlagG), ErrYmag);
    Res.MinAssymErr = MinErr;
    Res.AssymErr = LocalMinErr;
    Res.CooUnits = InPar.CooUnits;
        
    if (InPar.Plot)
        [Par,SPar] = polyfit(RefMag,Resid,InPar.PolyMagDeg);
        Xmag = (min(RefMag):InPar.StepMag:max(RefMag));
        Y   = polyval(Par,Xmag);
        semilogy(RefMag,Resid,'.');
        hold on
        semilogy(Xmag,Y,'k-')
        ErrY = SPar.normr^2/SPar.df;
        semilogy(Xmag,Y+ErrY.*5,'k--')
    end
    
    % Renormalize best fit parameters
    ParX = ParX./InPar.NormXY;
    ParY = ParY./InPar.NormXY;
    
    Res.TranC = populate_par(Res.TranC,Ix, ParX, ParErrX);
    Res.TranC = populate_par(Res.TranC,Iy, ParY, ParErrY);

    Res.ParX = ParX;
    Res.ParY = ParY;
    Res.ParErrX = ParErrX;
    Res.ParErrY = ParErrY;
    Res.NormXY  = InPar.NormXY;
    %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
    
    %Add the used flag to the residuals structure in order to save the used
    %objects.
    Res.FlagMag=Flag; 
    
    %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!    
    if all(InPar.BlockSize>=InPar.ImSize)
        % no Block statistics
    else
        
        if (~isempty(InPar.ImSize) && ~isempty(InPar.BlockSize))
            [ListEdge,ListCenter]=ImUtil.Im.image_blocks(InPar.ImSize,InPar.BlockSize,InPar.BufferSize,'simple');
            %NBy = floor(InPar.ImSize(2)./InPar.BlockSize(2));
            %NBx = floor(InPar.ImSize(1)./InPar.BlockSize(1));
            NBy = numel(unique(ListCenter(:,2)));
            NBx = numel(unique(ListCenter(:,1)));

            ListEdge = ListEdge - [InPar.ImSize(1).*0.5.*[1 1] InPar.ImSize(2).*0.5.*[1 1]];
            % convert ListEdge from pix to CooUnits
            
            ListEdge = ListEdge.*InPar.PixScale.*convert.angular('arcsec',InPar.CooUnits);
            Nblock = size(ListEdge,1);
            for Iblock=1:1:Nblock
                %
                Edge  = ListEdge(Iblock,:);
                
                FlagE = FlagG & (CatX>Edge(1) & CatX<Edge(2) & CatY>Edge(3) & CatY<Edge(4));

                Res.Block_rmsN(Iblock)  = Util.stat.rstd(Resid(FlagE));
                Res.Block_rmsXN(Iblock)  = Util.stat.rstd(ResidX(FlagE));
                Res.Block_rmsYN(Iblock)  = Util.stat.rstd(ResidY(FlagE));
                Res.Block_NsrcN(Iblock) = sum(FlagE);
            end
            Res.Block_rmsN  = reshape(Res.Block_rmsN,NBy,NBx);
            Res.Block_rmsXN = reshape(Res.Block_rmsXN,NBy,NBx);
            Res.Block_rmsYN = reshape(Res.Block_rmsYN,NBy,NBx);
            Res.Block_NsrcN = reshape(Res.Block_NsrcN,NBy,NBx);

    %         Res.worstBlock_rmsN = max(Res.Block_rmsN(:));
    %         Res.Block_N0        = sum(Res.Block_NsrcN(:)==0);
    %         Res.Block_meanN     = mean(Res.Block_NsrcN(:));
            %surface(Res.Block_rmsN)
        end
    end
end
    
    
