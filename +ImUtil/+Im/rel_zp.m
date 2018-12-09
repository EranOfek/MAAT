function Res=rel_zp(Mat,Err,varargin)
%--------------------------------------------------------------------------
% rel_zp function                                                   ImPhot
% Description: Calculate the best fit photometric zero point (ZP) per
%              image and mean magnitude per source to a set of instrumental
%              magnitudes.
%              The function receive a a matrix of the instrumental
%              magnitudes of the sources in a set of images. The matrix
%              size is Nsrc X Nepoch (i.e., row per source, column per
%              image). Using linear least squares, it solves for the best
%              fit ZP per image and mean magnitude per star.
%              I.e., it solves the equation M_ij = ZP_j + <M>_i, where M_ij
%              is the matrix of instrumental magnitudes, ZP_j is the ZP of
%              the j-th image, and <M>_i is the mean magnitude of the i-th
%              source. Given the calibrated magnitudes of some of the
%              sources the function  can calibrate the solution.
%              Furthermore, de-trending using additional parameters is
%              possible.
%              NOTE: This function is working but it is being improved.
% Input  : - Matrix of magnitudes of the sources in all the epochs.
%            The matrix size is [Nsrc X Nepoch].
%            Alternatively, this can be a structure array that contains the
%            matrix of instrumental magnitude (field name via the 'MagName'
%            parameter) as well as additional matrices like the error
%            matrix and any other matrices needed for de0-trending.
%          - Matrix of magnitude errors. If empty, then will look for the
%            errors in the first argument structure. If Errors are not
%            provided at all, the a default value will be used (see the
%            'DefErr' argument).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MagName'
%            'MagErrName'
%            'DefErr'
%            'JD'
%            'Sparse'
%            'Class'
%            'AddCol'
%            'CalibMag'
%            'CalibMagErr'
%            'DefCalibMagErr'
%            'AlgoLSQ'
%            'AlgoLSCOV'
%            'ConjGradPar'
%            'MagBinSize'         = 0.5;
%            'MinNobsBinMag'      = 3;
%            'MaxMedBinMag'       = 1;
%            'MaxStdBinMag'       = 0.3;
%            'InterpErrMethod'    = 'nearest';
%            'Niter'              = 2;
% Output : - A structure containing the best fit photometric zero points
%            and mean magnitudes.
%            The following fields are available:
%            .Par      - Vector of all fitted parameters.
%            .ParErr   - Vector of all fitted parameter errors.
%            .ResidAll -
%            .Chi2All  -
%            .Resid    -
%            .Chi2     -
%            .Ncal     -
%            .Npar     -
%            .Neq      -
%            .Nobs     -
%            .Ndof     -
%            .StdZP    -
%            .StdMeanMag-
%            .ZP       - Vector of relative zero points
%                        (Mag-Mag_ref).
%            .MeanMag  -
%            .AddPar   -
% Reference: Ofek et al. 2011 (Appendix A).
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Nep = 100; Nsrc=1000; Err=0.01; Mag=randn(Nsrc,Nep).*Err;
%          Err = Err.*ones(size(Mag));
%          ZP  = rand(1,Nep).*3;
%          MM  = rand(Nsrc,1).*7;
%          Mag = bsxfun(@plus,Mag,ZP);
%          Mag = bsxfun(@plus,Mag,MM);
%          Res = rel_zp(Mag,Err);
% Reliable: 
%--------------------------------------------------------------------------
import Util.fit.*
import timeseries.*

if (nargin==1)
    Err = [];
end


DefV.MagName            = {'MAG_PSF','MAG_APER_1','MAG_APER'};
DefV.MagErrName         = {'MAGERR_PSF','MAGERR_APER_1','MAGERR_APER'};
DefV.DefErr             = 0.05;
DefV.JD                 = [];   % if string then this is a field name in Mat
DefV.Sparse             = true;
DefV.Class              = 'double';
DefV.AddCol             = {};  %{'(mod(XWIN_IMAGE,1)-0.5).^2','(mod(YWIN_IMAGE,1)-0.5).^2','AIRMASS','COLOR'};
DefV.CalibMag           = [];
DefV.CalibMagErr        = [];
DefV.DefCalibMagErr     = 0.01;
DefV.AlgoLSQ            = 'lscov';   % '\','lscov','pcg','cgs','bicg',... (see ls_conjgrad.m)
DefV.AlgoLSCOV          = 'chol';    % 'chol'|'orth'
DefV.ConjGradPar        = {};
DefV.MagBinSize         = 0.5;
DefV.MinNobsBinMag      = 3;
DefV.MaxMedBinMag       = 1;
DefV.MaxStdBinMag       = 0.2;
DefV.InterpErrMethod    = 'nearest';
DefV.Niter              = 2;
DefV.VerbosePlot        = false;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.MagName))
    InPar.MagName = {InPar.MagName};
end
if (~iscell(InPar.MagErrName))
    InPar.MagErrName = {InPar.MagErrName};
end

% If Mat and Err are matrices then copy their content
% into a structure named Mat
if (isnumeric(Mat))
    InPar.MagName       = InPar.MagName{1};
    Tmp = Mat;
    clear Mat;
    Mat.(InPar.MagName) = Tmp;
else
    % assume Mat is a structure array - identify Magntude matrix in Mat
    FN = fieldnames(Mat);
    Nmagname = numel(InPar.MagName);
    IndName  = [];
    for Imagname=1:1:Nmagname
        if isempty(IndName) && any(strcmpi(InPar.MagName{Imagname},FN))
            IndName = Imagname;
        end
    end
    InPar.MagName = InPar.MagName{IndName};
end
% Size of mag matrix
% P - Number of epochs
% Q - Number of sources
[Nsrc,Nepoch] = size(Mat.(InPar.MagName));

if (~isempty(Err))
    % Assume Err is a numeric array
    InPar.MagErrName       = InPar.MagErrName{1};
    Mat.(InPar.MagErrName) = Err;
else
    % assume Mat contains Err
    FN = fieldnames(Mat);
    Nmagerrname = numel(InPar.MagErrName);
    IndName  = [];
    for Imagerrname=1:1:Nmagerrname
        if isempty(IndName) && any(strcmpi(InPar.MagErrName{Imagerrname},FN))
            IndName = Imagerrname;
        end
    end
    if (isempty(IndName))
        % Error wasn't found - set error to default value
        InPar.MagErrName = InPar.MagErrName{1};
        Mat.(InPar.MagErrName) = InPar.DefErr;
    else
        InPar.MagErrName = InPar.MagErrName{IndName};
    end
end


%--------------
%--- Get JD ---
%--------------
if (isempty(InPar.JD))
    JD = (1:1:Nepoch).';
else
    if (ischar(InPar.JD))
        JD = Mat.(InPar.JD);
    else
        JD = InPar.JD;
    end
end


%-------------------------
%--- The design matrix ---
%-------------------------
% construct the basic design matrix for the photometric
% calibration problem
if (isempty(InPar.CalibMag))
    [H]        = design_matrix_calib(Nepoch,Nsrc,InPar.Sparse,InPar.Class);
else
    [H,CalibH] = design_matrix_calib(Nepoch,Nsrc,InPar.Sparse,InPar.Class);
    H = [H;CalibH];
end

% The corresponding vector of observables contains
% a concatenation the vectors of sources of each image.
% Note that the matrix dimensions are [Nsrc X Nepoch]


Nes = Nepoch.*Nsrc;

% Adding the calibration block
if (isempty(InPar.CalibMag))
    % The following command generate a vector out of the matrix
    Y    = Mat.(InPar.MagName)(:);
    ErrY = Mat.(InPar.MagErrName)(:);
    
else
    Y             = zeros(Nepoch.*Nsrc + Nsrc,1);
    ErrY          = zeros(Nepoch.*Nsrc + Nsrc,1);
    Y(1:Nes)      = Mat.(InPar.MagName)(:);
    ErrY(1:Nes)   = Mat.(InPar.MagErrName)(:);
    % Calibration block was requested
    % Attach InPar.CalibMag to Y
    Y(Nes+1:end)  = InPar.CalibMag(:);
    
    % Attached InPar.CalibMagErr to ErrY
    if (isempty(InPar.CalibMagErr))
        % Set CalibMagErr to default value
        InPar.CalibMagErr = InPar.DefCalibMagErr.*ones(size(InPar.CalibMag));
    end
    ErrY(Nes+1:end) = InPar.CalibMagErr(:);
end


% Adding additional columns to the design matrix
if (~isempty(InPar.AddCol))
    error('AddCol option is not supported yet');
end


%----------------------------------
%--- Cleaning the design matrix ---
%----------------------------------
for Iiter=1:1:InPar.Niter
    % removing NaNs from the design matrix
    ErrY = ErrY.*ones(size(Y));
    Flag = all(~isnan(H),2) & ~isnan(Y) & ~isnan(ErrY);
    % FlagObs is like Flag but in which the calibration block is set to false
    FlagObs = Flag;
    FlagObs(Nes+1:end) = false;
    % H    = H(Flag,:);
    % Y    = Y(Flag);
    % ErrY = ErrY(Flag);

    Ncal = sum(Flag(Nes+1:end));
    Npar = size(H,2);
    Neq  = sum(Flag) - Ncal;
    Nobs = Neq - Ncal;
    Ndof = Nobs - Npar;


    %------------------------------------
    %--- Linear least square solution ---
    %------------------------------------
    switch lower(InPar.AlgoLSQ)
        case '\'
            Res.Par      = H(Flag,:)\Y(Flag);
            Res.ParErr   = nan(Npar,1);
            Res.ResidAll = Y - H*Res.Par;
            Res.Chi2All  = nan(size(Res.ResidAll));
            % The residuals are calculated only for the observations
            Res.Resid  = Y(FlagObs) - H(FlagObs,:)*Res.Par;
            Res.Chi2   = nan(size(Res.Resid));
        case 'lscov'
            [Res.Par,Res.ParErr] = lscov(H(Flag,:),Y(Flag),1./(ErrY(Flag).^2),InPar.AlgoLSCOV);
            Res.ResidAll = Y - H*Res.Par;
            Res.Chi2All  = sum((Res.ResidAll./ErrY).^2);
            % The residuals are calculated only for the observations
            Res.Resid  = Y(FlagObs) - H(FlagObs,:)*Res.Par;
            Res.Chi2   = sum((Res.Resid./ErrY(FlagObs)).^2);
        otherwise
            Res = ls_conjgrad(H,Y,ErrY,Inpar.AlgoLSQ,InPar.ConjGradPar{:});
            warning('Need to fix Res from ls_conjgrad.m')
    end

    % Add meta information
    if (~isempty(InPar.CalibMag))
        Res.ResidCalib = Y(Nobs+1:end) - H(Nobs+1:end,:)*Res.Par;
    end
    Res.Ncal  = Ncal;
    Res.Npar  = Npar;
    Res.Neq   = Neq;
    Res.Nobs  = Nobs;
    Res.Ndof  = Ndof;

    % Calculate the scatter per image ZP
    ResidMat  = reshape(Res.ResidAll,Nsrc,Nepoch);
    Res.StdZP = nanstd(ResidMat,[],1);
    Res.StdMeanMag = nanstd(ResidMat,[],2);



    % Calculate the scatter per source


    % Seperate the parameters
    Res.ZP      = Res.Par(1:Nepoch);
    Res.MeanMag = Res.Par(Nepoch+1:Nepoch+Nsrc);
    Res.AddPar  = Res.Par(Nepoch+Nsrc+1:end);


    OutCol           = {'MidBin','StartBin','EndBin',@numel,@median,@Util.stat.rstd};
%     [BinStd,OutCol]  = timeseries.binning([Res.MeanMag,Res.StdMeanMag],InPar.MagBinSize,[NaN NaN],OutCol,'astcat');
%     BinStdArith      = sprintf('numel>%d & median<%f & Util.stat.rstd<%f',InPar.MinNobsBinMag,InPar.MaxMedBinMag,InPar.MaxStdBinMag);
%     FlagBin          = logical(col_arith(BinStd,BinStdArith,'mat'));
%     MidBin           = col_get(BinStd,'MidBin');
%     MedBin           = col_get(BinStd,'median');
%     MedBin(~FlagBin) = NaN;

    [BinStd,OutCol]  = timeseries.binning([Res.MeanMag,Res.StdMeanMag],InPar.MagBinSize,[NaN NaN],OutCol,'mat');
    FlagBin = BinStd(:,4)>InPar.MinNobsBinMag & ...
              BinStd(:,5)<InPar.MaxMedBinMag & ...
              BinStd(:,6)<InPar.MaxStdBinMag;
    MidBin = BinStd(:,1);
    MedBin = BinStd(:,5);
    MedBin(~FlagBin) = NaN;


    % the next addition doesn't make sense:

    % ErrY             = ErrY+repmat(interp1(MidBin,MedBin,Res.MeanMag,InPar.InterpErrMethod),Nepoch,1);
    % ErrY             = [ErrY;InPar.CalibMagErr(:)];

    Mat.(InPar.MagErrName) = sqrt(bsxfun(@plus,Res.StdMeanMag.^2,Res.StdZP.^2) + Mat.(InPar.MagErrName).^2);
    ErrY = Mat.(InPar.MagErrName)(:);
    ErrY = [ErrY; InPar.CalibMagErr(:)];


    if (InPar.VerbosePlot)
        semilogy(Res.MeanMag,Res.StdMeanMag,'.')
        hold on;
    end
    
end
