function Fit=spec_fit1d_psf(Sim,SimSky,varargin)
%-------------------------------------------------------------------------
% spec_extract1 function                                           ImSpec
% Description: Given a traced and sky subtracted spectrum, extract the
%              intensity of a single spectrum, along the dispersion axis,
%              by fitting a 1-dimensional numerical PSF. The PSF is
%              constructed by binning the spectrum.
% Input  : - An image containing the sky subtracted spectrum in which
%            the spectrum dispersion direction is aligned with the X-axis.
%            This is a single image in one of the following forms:
%            (1) A structure array (SIM).
%            (2) A file name.
%            (3) A matrix.
%            (4) A cell array with a single file name string.
%            (5) A cell array with a single matrix image.
%          - A sky background image the same size of the first input
%            image. If this parameter is empty, then will attempt to
%            retrieve the background image from the 'BackImField' in
%            the first input argument.
%            Image types can be as the first input argument.
%            Default is empty.
%          * Arbitrary number of pairs of input arguments:
%            ...,keyword,value,..., where the keywords can be one of
%            the followings:
%            'DispDir'    - Dispersion direction of image {'x'|'y'}.
%                           Default is 'x'. 'y' option may not work in
%                           some cases.
%            'Center'   - Fixed position in the spatial direction in
%                         which the spectrum maximum is located.
%                         Default is the centeral pixel in the spatial
%                         position. If empty use default.
%            'BinDisp'  - Fit the intensity in bins along the dispersion
%                         axis, where this parameter set the semi-width
%                         in pixels in the bin.
%                         Default is 0 (i.e., bin size is 1)..
%            'BinDispPSF'- Construct the PSF in bins along the dispersion
%                         axis. If empty, then will use the entire spectrum
%                         to construct the PSF. If two element vector,
%                         then will use pixels in the specified dispersion
%                         direction range (i.e., [min max]).
%                         If this is a scalar then this is the semi-width
%                         around the current spectral point, in the
%                         dispersion direction to collapse in order to
%                         generate the PSF. Default is empty.
%            'CollapsePSF'- Method for collapsing data in bin
%                         {'median'|'mean'}. Default is 'median'.
%                         This is used for both the spectral binning
%                         and PSF collapsing.
%            'AperRad'  - Aperture radiuof PSF and aperture photometry,
%                         around trace center. Default is 6.
%                         If two elements are provided than [- +] radius
%                         of aperture.
%            'Gain'     - CCD gain for noise estimatation. Default is 1.
%            'RN'       - CCD readout noise [e-] for noise estimatation.
%                         Default is 10.
%            'Mask'     - Image bit mask to use in excluding bad pixels
%                         in the fitting process.
%                         If empty will attempt to use the 'Mask' field
%                         in the first input argument.
%            'MaskType' - If a mask doesn't exist then generate an empty
%                         mask of this class. Default is 'uint16'.
%            'RemoveBits'- List of indices of Mask bits which should not
%                         be used in the PSF fit. Default is [4 9].
%                         See def_bitmask_specpipeline.m
%            'MaxIter'  - Maximum number of rejection iterations.
%                         Default is 1. If 1 no rejection.
%            'Reject'   - Number of [Low, High] pixels to reject in
%                         each fit. Default is [0 1].
% Output : - Structure of fitted parameters:
%            .Disp  - Vector of indices for the dispersion axis.
%            .H     - Best fitted PSF heighted in units of the
%                     normalized PSF.
%            .ErrH  - Error in best fit PSF height.
%            .Mask  - Mask bit within aperture radius propogated to
%                     the 1D spectrum.
%            .SumPSF- The PSF normalization in each dispersion pixel.
%                     Negative normalizations should be ignored.
%            .Sum   - The integral of the best fit PSF.
%            .ErrSum- The error in the integral of the best fit PSF.
%            .RMS   - RMS of the best fit.
%            .Chi2  - \chi^2 of the best fit.
%            .Dof   - Number of degrees of freedom.
%            .AperSum- Aperture photometry sum.
%            .AperErr- Error in aperture photometry sum.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jan 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Fit=spec_fit1d_psf(BackSubSpec,[]);
% Reliable: 2
%-------------------------------------------------------------------------
import Util.array.*

ImageField  = 'Im';
%HeaderField = 'Header';
%FileField   = 'ImageFileName';
MaskField   = 'Mask';
BackImField = 'BackIm';
%ErrImField  = 'ErrIm';


if (nargin==1)
    SimSky = [];
end

% set default parameters
DefV.DispDir     = 'x';
DefV.Center      = [];
DefV.BinDisp     = 0;
DefV.BinDispPSF  = []; %300;   % if empty use entire image; if pair use range
DefV.CollapsePSF = 'median';
DefV.AperRad     = 6;
DefV.NormRad     = 6;
% noise
DefV.Gain        = 1;
DefV.RN          = 10;
% mask
DefV.Mask        = [];  % override the Mask field if available
DefV.MaskType    = 'uint16';
DefV.RemoveBits  = [4 9];
% fit
DefV.MaxIter     = 1;
DefV.Reject      = [0 1];

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

if (numel(InPar.AperRad)==1)
    InPar.AperRad = [InPar.AperRad, InPar.AperRad];
end
%if (numel(InPar.NormRad)==1),
%    InPar.NormRad = [InPar.NormRad, InPar.NormRad];
%end

% read image
Sim      = image2sim(Sim,varargin{:});
Image    = Sim.(ImageField);

% read background image (for noise calculation)
if (isempty(SimSky))
    % try to read background image from Sim
    if (isfield(Sim,BackImField))
        SkyImage = Sim.(BackImField);
    else
        % no background image - set to zero
        SkyImage = zeros(size(Image));
    end
else
    SimSky   = image2sim(SimSky,varargin{:});
    SkyImage = SimSky.(ImageField);     % read background from .Im field!
end

% read Mask image
if (~isempty(InPar.Mask))
    % InPar.Mask exist - override Sim.(MaskField)
    InPar.Mask = InPar.Mask;
else
    if (isfield(Sim,MaskField))
        InPar.Mask = Sim.(MaskField);
    else
        % set Mask to zeros
        InPar.Mask = zeros(size(Image),InPar.MaskType);
    end
end
 
% make sure dispersion axis is along the x-axis
switch lower(InPar.DispDir)
 case 'x'
    % do nothing
 case 'y'
    % rotate
    Image      = Image.';
    SkyImage   = SkyImage.';
    InPar.Mask = InPar.Mask.';
 otherwise
    error('Unknown DispAxis option');
end


% Spectrum size
[Nspat,Ndisp] = size(Image);
SpatCenter    = (Nspat+1).*0.5;   % centeral position of spatial axis
Fit.Disp      = (1:1:Ndisp).';    % vector of indices along the dispersion axis

% set Center to default
if (isempty(InPar.Center)),
   InPar.Center = SpatCenter;
end


AperInd = ((InPar.Center-InPar.AperRad(1)):1:(InPar.Center+InPar.AperRad(2)));
%NormInd = ((InPar.Center-InPar.NormRad(1)):1:(InPar.Center+InPar.NormRad(2)));

% set global psf
if (isempty(InPar.BinDispPSF))
    % set PSF range to entire image
    InPar.BinDispPSF = [1 Ndisp];
end
if (numel(InPar.BinDispPSF)==2)
    % assume BinDispPSF contains range
    % build PSF within range
    switch lower(InPar.CollapsePSF)
        case 'median'
            PSF = nanmedian(Image(AperInd,InPar.BinDispPSF(1):InPar.BinDispPSF(2)),2);
        case 'mean'
            PSF = nanmean(Image(AperInd,InPar.BinDispPSF(1):InPar.BinDispPSF(2)),2);
        otherwise
            error('Unknown CollapsePSF option');
    end
    % normalize PSF
    SumPSF = nansum(PSF);
    PSF = PSF./SumPSF;
    SumPSF = SumPSF.*ones(1,Ndisp);
end

%-------------------------
%--- Fit numerical PSF ---
%-------------------------
Fit.H         = zeros(Ndisp,1).*NaN;
Fit.ErrH      = zeros(Ndisp,1).*NaN;
Fit.SumPSF    = zeros(Ndisp,1).*NaN;
Fit.Sum       = zeros(Ndisp,1).*NaN;
Fit.ErrSum    = zeros(Ndisp,1).*NaN;
Fit.RMS       = zeros(Ndisp,1).*NaN;
Fit.Chi2      = zeros(Ndisp,1).*NaN;
Fit.Dof       = zeros(Ndisp,1).*NaN;
Fit.AperSum   = zeros(Ndisp,1).*NaN;
Fit.AperErr   = zeros(Ndisp,1).*NaN;

% the PSF is extracted from binning the spectrum
% with bin size BinDispN along the dispersion axis:
for Idisp=1:1:Ndisp,
   %-------------------
   %--- prepare PSF ---
   %-------------------
   
   if (numel(InPar.BinDispPSF)==1),
       % make a slidining window PSF
       % remove out of bound indices
       [MinInd,MaxInd] = Util.array.check_range(Ndisp,Idisp-InPar.BinDispPSF,Idisp+InPar.BinDispPSF);
   
       switch lower(InPar.CollapsePSF)
        case 'median'
           PSF = nanmedian(Image(AperInd,MinInd:MaxInd),2);
        case 'mean'
           PSF = nanmean(Image(AperInd,MinInd:MaxInd),2);
        otherwise
           error('Unknown CollapsePSF option');
       end
       % normalize PSF
       SumPSF(Idisp) = nansum(PSF);
       PSF           = PSF./SumPSF(Idisp);
   end
   
   %---------------------------------
   %--- Fit spectrum cut with PSF ---
   %---------------------------------

   %--------------------
   %--- Cut spectrum ---
   %--------------------
   % remove out of bound indices
   [MinInd,MaxInd] = Util.array.check_range(Ndisp,Idisp-InPar.BinDisp,Idisp+InPar.BinDisp);
   
   switch lower(InPar.CollapsePSF)
    case 'median'
       SpectrumCut = nanmedian(Image(AperInd,MinInd:MaxInd),2);
       SkyCut      = nanmedian(SkyImage(AperInd,MinInd:MaxInd),2);
    case 'mean'
       SpectrumCut = nanmean(Image(AperInd,MinInd:MaxInd),2);
       SkyCut      = nanmean(SkyImage(AperInd,MinInd:MaxInd),2);
    otherwise
       error('Unknown CollapsePSF option');
   end
   % create Mask cut
   MaskCut = sum_bitor(InPar.Mask(AperInd,MinInd:MaxInd),2);

   % Flag pixels that should be used in fit
   FlagUse = ~maskflag_check(MaskCut,sum(2.^(InPar.RemoveBits-1)),'and');
   
   % aperture photometry
   
   Fit.AperSum(Idisp) = sum(nanmedian(Image(AperInd,MinInd:MaxInd),2));
   AperBack           = sum(nanmedian(SkyImage(AperInd,MinInd:MaxInd),2));
   Fit.AperErr(Idisp) = (Fit.AperSum(Idisp).*InPar.Gain + AperBack.*InPar.Gain + InPar.RN.^2)./InPar.Gain;
   if (length(find(FlagUse))<2 || any(isnan(SpectrumCut(FlagUse))) )
       % do nothing
       Par    = NaN;
       ParErr = NaN;
       RMS    = NaN;
       Chi2   = NaN;
       Dof    = NaN;
   else
       Iiter = 0;
       while (Iiter<InPar.MaxIter)
          Iiter = Iiter + 1;
          % ignore PSF noise in fit:
          % The design matrix is 'PSF'
          % work in DN (rather than in e-)
          Y            = SpectrumCut(FlagUse);
          ErrY         = sqrt(abs(SpectrumCut(FlagUse)).*InPar.Gain + abs(SkyCut(FlagUse)).*InPar.Gain + InPar.RN.^2)./InPar.Gain;
          H            = PSF(FlagUse);
         
          [Par,ParErr] = lscov(H,Y,1./ErrY.^2);
          Resid        = Y - H*Par;
          RMS          = std(Resid);
          Chi2         = sum((Resid./ErrY).^2);
          Dof          = length(find(FlagUse)) - 1;  % 1 free par in fit

          % for next iteration - remove additional pixels
          
          if (sum(InPar.Reject)>0),
              [~, SortInd] = sort(Resid);
              if (InPar.Reject(1)>0),
                  FlagUse(SortInd(1:InPar.Reject(1))) = false;
              end
              if (InPar.Reject(2)>0),
                  FlagUse(SortInd(end-InPar.Reject(2)+1:end)) = false;
              end
          end
       end
   end
   Fit.H(Idisp)      = Par;
   Fit.ErrH(Idisp)   = ParErr;
   Fit.SumPSF        = SumPSF;
   Fit.Sum(Idisp)    = SumPSF(Idisp).*Fit.H(Idisp);
   Fit.ErrSum(Idisp) = SumPSF(Idisp).*Fit.ErrH(Idisp);
   Fit.RMS(Idisp)    = RMS;
   Fit.Chi2(Idisp)   = Chi2;
   Fit.Dof(Idisp)    = Dof;
end


%--- propogate mask ---
Fit.Mask = sum_bitor(InPar.Mask(AperInd,:),1);

