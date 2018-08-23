function [MatW,MatX,MatY,WaveCalibOut,RMS]=wavecalib_2d(ImageName,varargin)
%-----------------------------------------------------------------------
% waveclib_i function                                            ImSpec
% Description: Perform wavelength calibration for 2-D spectral images.
% Input  : - Image for which to preform Wavelengh calibration.
%            The dispersion should be along the x-axis.
%            If string, assume it to be a FITS file name.
%            If empty (e.g., []), prompt the user for FITS image name,
%            default is [].
%          *
%            'StartY'          - Spatial coordinate (along y-axis),
%                                in which to start the calibration.
%                                If empty (e.g., []), then the user will
%                                be asked to select this parameter
%                                interactively, default is [].
%                                If vector than perform the wavelength
%                                calibration along the specified spatial
%                                positions in this vector.
%            'WaveCalib'       - Wavelength calibration matrix, containing
%                                [dispersion_position(X), Wavelength(A)],
%                                default is no calibration matrix (i.e., [];
%                                In this case the wavelength calibration
%                                will be preform manually).
%            'CalibArcSpec'    - Calibrated Arc spectrum 
%                                [Wavelength(Ang), position(pixel), Intensity].
%                                default is [], if given then use
%                                wavecalib_1d_xcorr.m instead of
%                                wavecalib_1d_auto.m
%            'ArcType'         - Arc lines type (see get_arclines.m for more
%                                details) for automatic identification.
%                                Options are: {'Ar'|'Cd'|'Hg'|'Ne'|'Zn'|
%                                              'SA'|'SS'}, default is 'SS'.
%            'SigmaClip'       - Sigma clipping
%                                [-Number_sigma_low, +Number_sigma_high],
%                                for rejection of mis-identified lines,
%                                default is [-2, 2].
%            'PeakSearch'      - Semiwidth of line peak region size,
%                                default is 10 pixels.
%            'CenterMethod'    - Line peak centering method 
%                                {'max'|'wmean'|'stirling4'},
%                                default is 'stirling4'.
%            'WaveOrder'       - Polynomial order for fit in wavelength
%                                direction, default is 3.
%            'SpatialOrder'    - Polynomial order to fit in the spatial
%                                direction, default is 2.
%            'BinSizeSpatial'  - Bin size in spatial direction (y-axis),
%                                default is 3 (will be rounded to the
%                                nearest odd number).
%            'PlotLines'       - Plot identified arc lines {'y'|'n'},
%                                default is 'y'.
%            'PlotFit'         - Plot best fit for wavelength calibration,
%                                {'y'|'n'}, default is 'y'.
%            'InterpMethod'    - Interpolation method, default is 'linear'.
% Output : - Matrix of wavelength for each pixel in image.
%          - Matrix of X position for each pixel in image.
%          - Matrix of Y position for each pixel in image.
%          - Cell array of calibration for each pixel [X, Wavelength].
%            Each cell for each spatial position in which wavecalib
%            was attempted.
%          - Cell array of solutions RMS.
% Tested : Matlab 6.5
%     By : Eran O. Ofek            December 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------
if (nargin==0),
   ImageName = [];
end
Method = 'auto';     % default method for wavelengh calibration


% set default parameters
StartY           = [];
WaveCalib        = [];
CalibArcSpec     = [];
ArcType          = 'SS';
SigmaClip        = [-2, 2];
PeakSearch       = 10;
CenterMethod     = 'stirling4';
WaveOrder        = 3;
SpatialOrder     = 2;
BinSizeSpatial   = 3;
PlotLines        = 'y';
PlotFit          = 'y';
InterpMethod     = 'linear';

Narg = length(varargin);
for Iarg=1:2:Narg,
   switch lower(varargin{Iarg})
    case 'starty'
       StartY        = varargin{Iarg+1};
    case 'calibarcspec'
       CalibArcSpec  = varargin{Iarg+1};
    case 'wavecalib'
       WaveCalib     = varargin{Iarg+1};
    case 'arctype'
       ArcType       = varargin{Iarg+1};
    case 'sigmaclip'
       SigmaClip     = varargin{Iarg+1};
    case 'peaksearch'
       SPeakSearch   = varargin{Iarg+1};
    case 'centermethod'
       CenterMethod  = varargin{Iarg+1};
    case 'waveorder'
       WaveOrder     = varargin{Iarg+1};
    case 'spatialorder'
       SpatialOrder  = varargin{Iarg+1};
    case 'binsizespatial'
       BinSizeSpatial = varargin{Iarg+1};
    case 'plotlines'
       PlotLines      = varargin{Iarg+1};
    case 'plotfit'
       PlotFit        = varargin{Iarg+1};
    case 'interpmethod'
       InterpMethod   = varargin{Iarg+1};
    otherwise
       error('Unknown keyword option');
   end
end

% round to the nearest odd number
BinSizeSpatial = 2.*floor(BinSizeSpatial.*0.5) + 1;

if (isempty(ImageName)==1),
   % prompt the user for FITS image name
   ImageName = input('Enter FITS image name: ','s');
end

if (isstr(ImageName)==1),
   % read FITS file
   ImageMatrix = fitsread(ImageName); 
else
   error('Image must be FITS file - other options not supported yet');
   ImageMatrix = ImageName;
   clear ImageName;
end

[Ny,Nx]     = size(ImageMatrix);
[MatX,MatY] = meshgrid([1:1:Nx],[1:1:Ny]);

if (isempty(StartY)==1),
   % select StartY interactively
   disp('--- Open ds9 ---')
   % open FITS file in ds9 display   
   ds9_disp(ImageName,1);
   disp('--- Mark spatial position on image in which to start wavelength calibration');
   [Xspat,StartY] = ds9_getcoo(1,'image');
end

HalfBinSizeSpatial = 0.5.*(BinSizeSpatial - 1); 
VecStartY = [StartY(1)-HalfBinSizeSpatial: 1 :StartY(1)+HalfBinSizeSpatial];

% intensity as a function of dispersion position (uncalibrated
% wavelength=pixels) at the StartY spatial position.
IntensityDispDir_AtY{1} = median(ImageMatrix(VecStartY,:),1).';
VecX                    = [1:1:Nx].';

IntensityDispDir_AtY{1}

% If preliminary line lists exist call: wavecalib_1d_auto.m
% If calibrated spectrm exist call: wavecalib_1d_xcorr.m
% else call: wavecalib_1d_int.m, when exit call: wavecalib_1d_auto.m

%---------------------------------------------------
%--- First wavelength calibration in 2d spectrum ---
%---------------------------------------------------
if (size(WaveCalib,1)>1),
   % Enough lines in initial WaveCalib: use wavecalib_1d_auto.m
   [SelWaveCalib,WaveCalibOut{1},RMS{1},ParWaveCalib{1}] = wavecalib_1d_auto(IntensityDispDir_AtY{1},...
                                 'Ignore','y',varargin{:});
   Method = 'auto';
elseif (isempty(CalibArcSpec)==0),
   % CalibArcSpec exits: use wavecalib_1d_xcorr.m
   error('xcorr option not implemented yet');
   Method = 'xcorr';
else
   % use interactive wavelength calibration
   [SelWaveCalib,WaveCalibOut{1},RMS{1},ParWaveCalib{1}] = wavecalib_1d_int(IntensityDispDir_AtY{1},...
                                 'Ignore','y',varargin{:});

end
%----------------------------------
%--- Define the solution matrix ---
%----------------------------------
Xw(1:Nx,1) = WaveCalibOut{1}(:,1);          % X pixel position
Yw(1:Nx,1) = StartY(1).*ones(size(VecX));   % Y pixel position
Ww(1:Nx,1) = WaveCalibOut{1}(:,2);          % corresponding wavelength


%---------------------------------------------------------
%--- run wavecalib in additional spatial positions (Y) ---
%---------------------------------------------------------
% for each spatial position
for Isp=2:1:length(StartY),

   VecStartY = [StartY(1)-HalfBinSizeSpatial: 1 :StartY(1)+HalfBinSizeSpatial];
   % intensity as a function of dispersion position (uncalibrated
   % wavelength=pixels) at the StartY spatial position.
   IntensityDispDir_AtY{Isp} = median(ImageMatrix(VecStartY,:),1).';
   %VecX                    = [1:1:Nx].';    % already defined

   switch lower(Method)
    case 'auto'
       % automatic wave calib
       [SelWaveCalib,WaveCalibOut{Isp},RMS{Isp},ParWaveCalib{Isp}] = wavecalib_1d_auto(IntensityDispDir_AtY{Isp},...
                                 'Ignore','y',varargin{:},...
                                 'WaveCalib',SelWaveCalib);

    case 'xcorr'
       error('xcorr option not implemented yet');

    otherwise
       error('Unknown Method option');
   end

   %----------------------------------
   %--- Add to the solution matrix ---
   %----------------------------------
   Xw(1:Nx,Isp) = WaveCalibOut{Isp}(:,1);          % X pixel position
   Yw(1:Nx,Isp) = StartY(Isp).*ones(size(VecX));   % Y pixel position
   Ww(1:Nx,Isp) = WaveCalibOut{Isp}(:,2);          % corresponding wavelength

end

% Output is a matrix of wavelength per each pixel position.
MatW = griddata(Xw,Yw,Ww,MatX,MatY,InterpMethod);
