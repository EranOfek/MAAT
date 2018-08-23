function FitRes=fit_transform(AstC,AstR,varargin)
%--------------------------------------------------------------------------
% fit_transform function                                     class/@AstCat
% Description: OBSOLETE: use ImUtil.pattern.fit_transform instead
% Input  : - An AstCat object.
%          - An AstCat object containing a single reference catalog.
%            The coordinates of this catalog are used as the "observed
%            position".
%            If empty, then will use one of the catalog in the first
%            input AstCat object as a reference catalog (see 'SelRef').
%            Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'AstT'   - 
%            'ColCatXY'- Cell array of column names containing the
%                         X and Y column names.
%                         Default is {'XWIN_IMAGE','YWIN_IMAGE'}.
%            'SelRef' - Method for selecting refrence catalog from the
%                       input catalogs. Options are:
%                       'first' - use first catalog in the first input
%                                 argument.
%                       'last'  - use last catalog in the first input
%                                 argument.
%                       'maxn'  - use the catalog with the maximum number
%                                 of sources. Default.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstR = AstCat; AstR.Cat=rand(1000,2).*1000;
%          AstR.ColCell = {'XWIN_IMAGE','YWIN_IMAGE'};
%          AstR = colcell2col(AstR);
%          AstC = AstR; AstC.Cat = AstR.Cat + 5 + randn(1000,2);
%          F=fit_transform(AstC,AstR,'ColRefXY',{'XWIN_IMAGE','YWIN_IMAGE'})
% Reliable: 2
%--------------------------------------------------------------------------
import timeseries.*

DefV.AstT               = [];
DefV.ColCatXY           = {'XWIN_IMAGE','YWIN_IMAGE';'X','Y'};  % If ColNames empty & ~AstCat -> [1 2]
DefV.ColRefXY           = {'XWIN_IMAGE','YWIN_IMAGE';'X','Y'}; % {'X','Y'};  
DefV.ColRefMag          = {'Mag','MAG_PSF','MAG_APER'};
DefV.CatSelect          = {};    %{'SN',10,Inf};
DefV.ColRefColor        = 'Color';
DefV.SelRef             = 'maxn';
DefV.Nsigma             = 3;      % robust sigma clipping % If inf - no clip
DefV.NormXY             = 1;
DefV.AlgoLSQ            = 'chol';  %'cgs'; %'pcg';  %'chol';
DefV.Tol                = 1e-12;
DefV.MagBin             = 0.5;
DefV.MinNinBin          = 3;
DefV.MagInterpMethod    = 'linear';
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% Select the ColCatXY, ColRefXY colum names:
[InPar.ColCatXY] = select_exist_colnames(AstC,InPar.ColCatXY);
[InPar.ColRefXY] = select_exist_colnames(AstR,InPar.ColRefXY);

% Select the ColRefMag column name:
[InPar.ColRefMag] = select_exist_colnames(AstR,InPar.ColRefMag(:));
% Convert to a string
InPar.ColRefMag   =  InPar.ColRefMag{1};

Ncat = numel(AstC);

if (nargin==1)
    AstR = [];
end

if (isempty(AstR))
    % If AstR is not provided then use one of the AstC
    % as a reference catalog.
    switch lower(InPar.SelRef)
        case 'first'
            RefInd = 1;
        case 'last'
            RefInd = Ncat;
        case 'maxn'
            Nsrc   = sizecat(AstC);
            [~,RefInd] = max(Nsrc);
        otherwise
            error('Unknown SelRef option');
    end
    AstR = AstC(RefInd);
end

% define the transformation
if (isempty(InPar.AstT))
    % default is shift transformation"
    %CellTran = {'x_shift',[];'y_shift',[];'x_rot',[];'y_rot',[]};
    CellTran = {'x_shift',[];'y_shift',[];'x_rot',[];'y_rot',[];'x_platetilt',[];'y_platetilt',[]};
    %CellTran = {'x_shift',[];'y_shift',[];'x_rot',[];'y_rot',[];'x_p2',[];'y_p2',[]};
    %CellTran = {'x_shift',[];'y_shift',[];'x_rot',[];'y_rot',[];'x_platetilt',[];'y_platetilt',[];'x_ofad',[];'y_ofad',[]};
    
    AstT     = cell2asttran(CellTran);
else
    if (isasttran(InPar.AstT))
        AstT = InPar.AstT;
    elseif (iscell(InPar.AstT))
        AstT = cell2asttran(InPar.AstT);
    else
        error('Illegal AstT option');
    end
end


FitRes = struct_def({'Par','ParErr','Resid','ResidX','ResidY','ResidT',...
                     'Rrms','Chi2','Nobs','Npar','Ndof','rms','AstT','Flag','rmsAll','RrmsAll','Wrms','MagRef','Xref','Yref'},...
                    Ncat,1);
for Icat=1:1:Ncat
    % for each input catalog
    
    % Select rows in catalog
    [AstC(Icat),FlagSelect] = row_select(AstC(Icat),InPar.CatSelect);
    
    % construct the design matrix for the fit:
    %[~,Design] = transform(AstT,AstC(Icat).Cat,'ColXY',InPar.ColCatXY,'NormXY',InPar.NormXY); %,...
    [~,Design] = transform(AstT,AstC(Icat),'ColXY',InPar.ColCatXY,'NormXY',InPar.NormXY); %,...
                                              % 'NormX',InPar.NormX,'NormY',InPar.NormY);
    % Flag is set to false for bad sources (NaN):
    FlagXY = Design.Flag;
    
    
    ColIndR = colname2ind(AstR,InPar.ColRefXY);
    Xref    = AstR.Cat(FlagSelect,ColIndR(1)); %./InPar.NormX;
    Yref    = AstR.Cat(FlagSelect,ColIndR(2)); %./InPar.NormY;
    
    XYref   = [Xref; Yref];
    Err     = ones(size(XYref));
    Nsrc    = numel(Xref);
    % Magnitude of sources in reference catalog
    % This is needed in order to calculate the rms as a function of
    % magnitude
    if (isempty(InPar.ColRefMag))
        MagRef = [];
    else
        MagRef  = col_get(AstR,InPar.ColRefMag);
        MagRef  = MagRef(FlagSelect);
    end
    
    % Fit - first iteration:
    
    %Norm = 1./max(abs(Design.H));
    %H    = bsxfun(@times,Design.H,Norm);
    
    switch lower(InPar.AlgoLSQ)
        case {'chol','orth'}
            [FitRes(Icat).Par,FitRes(Icat).ParErr] = lscov(Design.H(FlagXY,:),XYref(FlagXY),1./(Err(FlagXY).^2),InPar.AlgoLSQ);
        otherwise
            [FitRes(Icat).Par,FitRes(Icat).ParErr] = ls_conjgrad(Design.H(FlagXY,:),XYref(FlagXY),1./(Err(FlagXY).^2),InPar.AlgoLSQ,InPar.Tol);
    end
    
    %[FitRes(Icat).Par,FitRes(Icat).ParErr] = lscov(H,XYref,1./(Err.^2));
%     %FitRes(Icat).Par = FitRes(Icat).Par./Norm.';
    
    FitRes(Icat).Par          = FitRes(Icat).Par./Design.VecNorm.';
    FitRes(Icat).Resid        = XYref - Design.H*FitRes(Icat).Par;
    FitRes(Icat).ResidX       = FitRes(Icat).Resid(1:Nsrc);
    FitRes(Icat).ResidY       = FitRes(Icat).Resid(Nsrc+1:end);
    FitRes(Icat).ResidT       = sqrt(FitRes(Icat).ResidX.^2 + FitRes(Icat).ResidY.^2);
    FitRes(Icat).MagRef       = MagRef;
    FitRes(Icat).Xref         = XYref(1:Nsrc);
    FitRes(Icat).Yref         = XYref(Nsrc+1:end);
    Nx = numel(FitRes(Icat).ResidX);
    % Robust rms
    FitRes(Icat).Rrms         = rstd(FitRes(Icat).ResidT);
    
    % Estimate the std as a function of magnitude
    if (~isempty(MagRef))
        BinRMS = binning([MagRef,FitRes(Icat).ResidT],InPar.MagBin,[NaN NaN],{'MidBin', @numel, @median});
        BinRMS = BinRMS(BinRMS(:,2)>InPar.MinNinBin,:);
        Err    = interp1(BinRMS(:,1),BinRMS(:,3),MagRef,InPar.MagInterpMethod);
        Err    = repmat(Err,2,1);
    end
    
    % Fits - second iteration
    if (~isinf(InPar.Nsigma))
        Flag   = abs(FitRes(Icat).ResidT)<(FitRes(Icat).Rrms.*InPar.Nsigma);
        FlagXY = FlagXY & [Flag;Flag];
        FlagXY  = FlagXY & ~isnan(Err); %   [MagRef;MagRef]>7 & [MagRef;MagRef]<11;
        %Nsrc   = sum(Flag);
        switch lower(InPar.AlgoLSQ)
            case {'chol','orth'}
                [FitRes(Icat).Par,FitRes(Icat).ParErr] = lscov(Design.H(FlagXY,:),XYref(FlagXY),1./(Err(FlagXY).^2),InPar.AlgoLSQ);
            otherwise
                [FitRes(Icat).Par,FitRes(Icat).ParErr] = ls_conjgrad(Design.H(FlagXY,:),XYref(FlagXY),1./(Err(FlagXY).^2),InPar.AlgoLSQ,InPar.Tol);
        end   
        FitRes(Icat).Par          = FitRes(Icat).Par./Design.VecNorm.';
        FitRes(Icat).Resid        = XYref - Design.H*FitRes(Icat).Par;
        FitRes(Icat).ResidX       = FitRes(Icat).Resid(1:Nsrc);
        FitRes(Icat).ResidY       = FitRes(Icat).Resid(Nsrc+1:end);
        FitRes(Icat).ResidT       = sqrt(FitRes(Icat).ResidX.^2 + FitRes(Icat).ResidY.^2);
    else
        Flag   = true(Nx,1);
        FlagXY = true(2.*Nx,1);
    end
    Flag = FlagXY(1:Nsrc);
    FitRes(Icat).Chi2   = sum((FitRes(Icat).Resid(FlagXY)./Err(FlagXY)).^2);
   
    FitRes(Icat).Nobs   = sum(FlagXY); %numel(FitRes(Icat).Resid);
    FitRes(Icat).Npar   = numel(FitRes(Icat).Par);
    FitRes(Icat).Ndof   = FitRes(Icat).Nobs - FitRes(Icat).Npar;
    % regular rms (std)
    FitRes(Icat).rmsAll = nanstd(FitRes(Icat).ResidT);
    % robust rms
    FitRes(Icat).RrmsAll= rstd(FitRes(Icat).ResidT);
    FitRes(Icat).rms    = std(FitRes(Icat).ResidT(Flag));
    FitRes(Icat).Rrms   = rstd(FitRes(Icat).ResidT(Flag));
    % weighted rms
    ErrX = Err(1:Nsrc);
    [~,~,FitRes(Icat).Wrms]=wmean([FitRes(Icat).ResidT(Flag),ErrX(Flag)]);

    % populate best fit parameters in transformation
    FitRes(Icat).AstT   = populate_par(AstT,Design.ParIndex,FitRes(Icat).Par);
    FitRes(Icat).Flag   = Flag;
    % Some statistics
    
    % Calculate the number of ancor points in each section of the CCD
    
    
    % plots
    %semilogy(MagRef,FitRes(Icat).ResidT,'.')
end




