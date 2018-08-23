function [Fit]=spec_extract1(Image,SkyImage,varargin);
%-------------------------------------------------------------------------
% spec_extract1 function                                           ImSpec
% Description: Given a traced and sky subtracted spectrum, extract the
%              intensity of a single spectrum, along the dispersion axis.
% Input  : - A 2-dimensional matrix, or a fits file, containing a traced
%            and sky subtracted spectrum.
%          - A 2-dimensional matrix, or a fits file, containing the sky,
%            in the original sky subtrcated spectrum (See second
%            output of spec_skysub.m).
%            The sky is used fo the fitting schemes.
%            Default is empty matrix (i.e., []), if empty then set sky
%            value to zero.
%          * Arbitrary number of pairs of input arguments:
%            ...,keyword,value,..., where the keywords can be one of
%            the followings:
%            'Method'   - Fitting spectrum intensity method:
%                         'Gauss' - Fit a Gaussian with unknown height
%                                   and width.
%                         'GaussC'- Fit a Gaussian with unknown height
%                                   width, and center.
%                         'Num'   - Fit a numerical PSF, with unknown
%                                   hight (default).
%                         'Aper'  - Sum the intensity within an
%                                   aperture around the trace center.
%                                   (See 'Aper' keyword).
%                                   Note that if BinDisp>0, then the
%                                   data is first binned and only then
%                                   summed (i.e., to remove CR).
%            'Nrej'     - Number of [low, high] pixels to reject in
%                         each fit, default is [0 1].
%            'Center'   - Fixed position in the spatial direction in
%                         which the spectrum maximum is located.
%                         Default is the centeral pixel in the spatial
%                         position.
%            'BinDisp'  - Fit the intensity in bins along the dispersion
%                         axis, where this parameter set the semi-width
%                         in pixels in the bin,
%                         default is 0 (i.e., bin size is 1)..
%            'BinDispC' - In the 'GauusC' method, find the center of the
%                         Gaussian in bins along the dispersion axis
%                         (in order to increase S/N). This parameter set
%                         the semi-width in pixels of the bin,
%                         default is 20
%            'BinDispW' - In the 'Gauus' and 'GaussC' methods, find the
%                         width of the Gaussian in bins along the
%                         dispersion axis (in order to increase S/N).
%                         This parameter set the semi-width in pixels
%                         of the bin, default is 50.
%            'BinDispN' - In the 'Num' method, construct the PSF in bins
%                         along the dispersion axis.
%                         This parameter set the semi-width in pixels
%                         of the bin, default is 50.
%            'BinCol'   - Method for collapsing data in bin
%                         {'median'|'mean'}, default is 'median'.
%            'Aper'     - The aperture lowest and highest pixels
%                         [low, high] along the dispersion axis.
%                         In the 'Aper' method return the sum of the
%                         intensity within the aperture.
%                         In the other methods, return the sum of the
%                         intensity in the fitted function, within the
%                         aperture. Default aperture is central pixel +/-3.
%                         Note that this also indicate the fitted region
%                         (i.e., pixels outside the aperture are ignored in
%                         the fitting process).
%            'Gain'     - CCD gain used in the fitting S/N estimation.
%                         Default is 1.
%            'RN'       - CCD readout noise [electrons] used in the
%                         fitting S/N estimation. Default is 6.
%            'DispAxis' - Dispersion axis {'x'|'y'}, default is 'x'.
% Output : - Structure of fitted parameters -
%            For all methods:
%            .Disp - Vector of indices for the dispersion axis.
%            .Sum  - Sum of counts in aperture. Fot the functions
%                    fitting this is the integal of the function
%                    within the aperture [counts].
%            .ErrSum- Error on sum inside aperture [counts].
%            .SkyAper - Sum of sky within aperture [counts].
%            .SN   - Signal/Noise.
%            .H    - Best fit height. For the 'Aper' method, this
%                    is the value of the central pixel [counts].
%            .ErrH - Error in best fit Gaussian height.
%            .W    - Best fit Gaussian width
%            .ErrW - Error in best fit Gaussian width
%            .C    - Best fit Gaussian center
%            .ErrC - Error in best fit Gaussian center
%            'Chi2'  - Chi2 of best fit.
%            'Dof'   - Number of degrees of freedoms.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                   January 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-------------------------------------------------------------------------


if (isstr(Image)==1),
   % Spec is a fits file
   Spec = fitsread(Image);
else
   Spec = Image;
   clear Image;
end

if (isempty(SkyImage)==1),
   Sky = zeros(size(Spec));
else
   if (isstr(SkyImage)==1),
      % Spec is a fits file
      Sky = fitsread(Image);
   else
      Sky = SkyImage;
      clear SkyImage;
   end
end

% set default parameters
Method    = 'Num';
Nrej      = [0 1];
Center    = [];
BinDisp   = 0;
BinDispC  = 20;
BinDispW  = 50;
BinDispN  = 50;
BinCol    = 'median';
Aper      = [];
Gain      = 1;
RN        = 6;
DispAxis  = 'x';

Narg      = length(varargin);
for Iarg=1:2:Narg-1,
   switch lower(varargin{Iarg})
    case 'method'
       Method      = varargin{Iarg+1};
    case 'nrej'
       Nrej        = varargin{Iarg+1};
    case 'center'
       Center      = varargin{Iarg+1};
    case 'bindisp'
       BinDisp     = varargin{Iarg+1};
    case 'bindispc'
       BinDispC    = varargin{Iarg+1};
    case 'bindispw'
       BinDispW    = varargin{Iarg+1};
    case 'bindispn'
       BinDispN    = varargin{Iarg+1};
    case 'bincol'
       BinCol      = varargin{Iarg+1};
    case 'aper'
       Aper        = varargin{Iarg+1};
    case 'gain'
       Gain        = varargin{Iarg+1};
    case 'rn'
       RN          = varargin{Iarg+1};
    case 'dispaxis'
       DispAxis    = varargin{Iarg+1};
    otherwise
       error('Unknown keyword option');
   end
end

switch lower(DispAxis)
 case 'x'
    % do nothing
 case 'y'
    % rotate
    Spec   = Spec.';
 otherwise
    error('Unknown DispAxis option');
end

% Spectrum size
[Nspat,Ndisp] = size(Spec);
SpatCenter    = round((Nspat+1).*0.5);  % centeral position of spatial axis
Fit.Disp      = [1:1:Ndisp].';          % vector of indices along the dispersion axis

% set Center to default
if (isempty(Center)==1),
   Center = SpatCenter;
end

% set Aper to default
if (isempty(Aper)==1),
   Aper   = [Center - 3, Center + 3];
end

% Vector of indices in Spec (spatial direction) which are inside aperture:
AperVec = [Aper(1):1:Aper(2)].';

% Index of Center in the AperVec vector:
AperCenter = round((Aper(2)-Aper(1)).*0.5);

switch Method
 case 'Aper'
    if (BinDisp==0),
       %----------------------------------
       %--- simple aperture extraction ---
       %----------------------------------
       Fit.Sum       = sum(Spec(AperVec,:),1).';
       Fit.H         = Spec(Center,:).';
       Fit.SkyAper   = sum(Sky(AperVec,:),1).';
       Fit.SN        = Fit.Sum.*Gain./sqrt((Fit.Sum + Fit.SkyAper).*Gain + RN.^2);
       Fit.ErrSum    = Fit.Sum./Fit.SN;
    else
       % Bin the data before sum
       switch lower(BinCol)
        case 'median'
           % Median filtering along the dispersion axis
           % each pixel
           Nbin      = 1 + BinDisp.*2;
           MedFilt   = medfilt1(Spec(AperVec,:), Nbin, [], 2);
           Fit.H         = MedFilt(AperCenter,:).';
           Fit.Sum       = sum(MedFilt).';
           Fit.SkyAper   = sum(Sky(AperVec,:),1).';
           Fit.SN        = Nbin.*Fit.Sum.*Gain./sqrt((Nbin.*Fit.Sum + Fit.SkyAper).*Gain + RN.^2);
           Fit.ErrSum    = Fit.Sum./Fit.SN;
        case 'mean'
           error('Filter mean is not implemented');
        otherwise
           error('Unknown BinCol option');
       end
    end
 otherwise
    %------------------------------------
    %--- prepare spectrum cut for fit ---
    %------------------------------------

    switch Method
     case 'Num'
        %-------------------------
        %--- Fit numerical PSF ---
        %-------------------------
        % the PSF is extracted from binning the spectrum
        % with bin size BinDispN along the dispersion axis:

        FitPSF = spec_fit1d_psf(Spec,Sky,'Nrej',Nrej,'Center',Center,'BinDisp',BinDisp,'BinDispN',BinDispN,'BinCol',BinCol,'Aper',Aper,'Gain',Gain,'RN',RN,'DispAxis',DispAxis);

        Fit.H      = FitPSF.H;
        Fit.ErrH   = FitPSF.ErrH.';
        Fit.Sum    = FitPSF.Sum.';
        Fit.ErrSum = FitPSF.ErrSum.';
        Fit.Chi2   = FitPSF.Chi2.';
        Fit.Dof    = FitPSF.Dof.';

     case 'Gauss'
        error('Gauss not implemented');
     case 'GaussC'
        error('GaussC not implemented');
     otherwise
        error('Unknwon Method option');
    end
end
