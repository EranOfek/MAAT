function Res=calib_flux_lsq(Mag,Err,varargin)
% Best fit zeropoints and mean magnitudes to a matrix of magnitudes.
% Package: timeseries
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
% Input  : - Matrix of instrumental magnitudes [Src,Epoch] to calibrate.
%            If no arguments are provided than will run in simulation mode.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=timeseries.calib_flux_lsq;
% Reliable: 
%--------------------------------------------------------------------------

if (nargin<2)
    Err = [];
end

if (nargin==0)
    Nsrc = 100;
    Nep  = 300;
    
    ZP      = 3.*rand(1,Nep);
    MeanMag = rand(Nsrc,1).*5 + 15;
    Err     = 0.01;
    
    Mag     = MeanMag*ZP;
    Mag     = Mag + Err.*randn(size(Mag));
    
end
    

DefV.Niter                = 2;
DefV.FlagSrc              = true;
DefV.CalibMag             = nan(0,0);
DefV.CalibErr             = 1e-2;
DefV.rmsflux_selectPar    = {};
DefV.Sparse               = true;
DefV.Class                = 'double';
DefV.FitMethod            = 'lscov';
DefV.ConjGradMethod       = 'pcg';
DefV.lscov_Algo           = 'chol';    % 'chol' | 'orth'
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (isempty(Err))
    Err = 1;
end

[Nsrc,Nep] = size(Mag);
if (numel(Err)==1)
    Err = Err+zeros(Nsrc,Nep);
end


% Prepare the design matrix and observation vectors
[H,CalibH] = Util.fit.design_matrix_calib(Nsrc,Nep,InPar.Sparse,InPar.Class);
MagVec     = Mag(:);
ErrVec     = Err(:);

% Prepare initial flag of sources to use in the solution
FlagSrc = InPar.FlagSrc;
if (numel(FlagSrc)==1)
    FlagSrc = true(Nsrc,1);
else
    if (numel(FlagSrc)~=Nsrc)
        error('Number of elements in FlagSrc must be 1 or equal to the number of sources');
    end
end

% FlagSrc is per star, so need to
% duplicate FlagSrc for all images
FlagSrc = repmat(FlagSrc,Nep,1);

% By default there are no calibration stars
FlagSrcCalib = [];

% add calibration block into the design matrix
if (~isempty(InPar.CalibMag))
    % add calibration part into the design matrix
    if (numel(InPar.CalibMag)==Nsrc)
        FlagNN = ~isnan(InPar.CalibMag);
        Hall = [H;CalibH(FlagNN,:)];
        MagVec = [MagVec; InPar.CalibMag(FlagNN)];

        % Flag for calibration sources
        FlagSrcCalub         = false(Nsrc,1);
        FlagSrcCalib(FlagNN) = true;
        
        if (numel(InPar.CalibErr)==1)
            InPar.CalibErr = InPar.CalibErr.*ones(Nsrc,1);
            ErrVec = [ErrVec; InPar.CalibErr(FlagNN)];
        end
        
    else
        error('Number of calibration sources must be either 0 or equal to the number of sources');
    end
else
    % in case we are not using the calibration part set Hall=H
    Hall = H;
end

% fitting / iterations
for Iiter=1:1:InPar.Niter

    % clean list of sources to use in the fit
    if (Iiter>1)
        % after the first iteration select stars using various methods
        
        [Flag,ResS]=timeseries.rmsflux_select(Res.MeanFlux,Res.StarsRelRMS,InPar.rmsflux_selectPar{:}); %'Plot',true);
        
    else
        % in the first iteration select all stars
        FlagSrcIter = FlagSrc;
    end
    % Flag sources including the calibration part
    FlagSrcAll = [FlagSrcIter; FlagSrcCalib];
    
    switch lower(InPar.FitMethod)
        case 'lscov'
            [Par,ParErr] = lscov(Hall(FlagSrcAll,:),MagVec(FlagSrcAll),1./(ErrVec(FlagSrcAll).^2),InPar.lscov_Algo);
        case '\'
            Par = Hall(FlagSrcAll,:)\MagVec(FlagSrcAll);
        case 'conjgrad'
            [Par,ParErr] = Util.fit.ls_conjgrad(Hall(FlagSrcAll,:),MagVec(FlagSrcAll),ErrVec(FlagSrcAll),InPar.ConjGradMethod);
        otherwise
            error('Unknown FitMethod');
    end
    
    % note that the residuals are calculated without the calibration part
    Resid = MagVec - H*Par;
    RMS   = std(Resid);
    Chi2  = sum((Resid./ErrVec).^2);
    Neq   = sum(FlagSrc);       % number of equations
    Ncon  = sum(FlagSrcCalib);  % number of constraints
    Npar  = Nsrc + Nep;
    Ndof  = Neq - Npar;
end


Res.Par   = Par;
Res.Resid = Resid;
Res.RMS   = RMS;
Res.Chi2  = Chi2;
Res.Neq   = Neq;
Res.Ncon  = Ncon;
Res.Npar  = Npar;
Res.Ndof  = Ndof;

