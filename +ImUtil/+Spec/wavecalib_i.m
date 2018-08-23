function [SelWaveCalib,WaveCalibOut]=wavecalib_i(ImageName,varargin)
%-----------------------------------------------------------------------
% waveclib_i function                                            ImSpec
% Description: Perform wavelength calibration for spectral images.
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
%            'WaveCalib'       - Wavelength calibration matrix, containing
%                                [dispersion_position(X), Wavelength(A)],
%                                default is no calibration matrix (i.e., [];
%                                In this case the wavelength calibration
%                                will be preform manually).
%                                
%            'PeakSearch'      - Semiwidth of line peak region size,
%                                default is 3 pixels.
%            'CenterMethod'    - Line peak centering method {'max'|'wmean'},
%                                default is 'wmean'.
%            'WaveOrder'       - Polynomial order for fit in wavelength
%                                direction, default is 3.
%            'SpatialOrder'    - Polynomial order to fit in the spatial
%                                direction, default is 2.
%            'BinSizeSpatial'  - Bin size in spatial direction (y-axis),
%                                default is 3 (will be rounded to the
%                                nearest odd number).
%            'InterpMethod'    - Interpolation method, default is 'linear'.
% Output : - List of the peak wavelength for the selected line peaks
%            [X, Wavelength, ErrX, Y, FlagMax].
%          - Calibration matrix for each pixel [X, Wavelength].
% Tested : Matlab 6.5
%     By : Eran O. Ofek            December 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------
if (nargin==0),
   ImageName = [];
end

% set default parameters
StartY           = [];
WaveCalib        = [];
PeakSearch       = 3;
CenterMethod     = 'wmean';
WaveOrder        = 3;
SpatialOrder     = 2;
BinSizeSpatial   = 3;
InterpMethod     = 'linear';

Narg = length(varargin);
for Iarg=1:2:Narg,
   switch lower(varargin{Iarg})
    case 'starty'
       StartY        = varargin{Iarg+1};
    case 'wavecalib'
       WaveCalib     = varargin{Iarg+1};
    case 'peaksearch'
       SPeakSearch   = varargin{Iarg+1};
    case 'centermethod'
       CenterMethod  = varargin{Iarg+1};
    case 'waveorder'
       WaveOrder      = varargin{Iarg+1};
    case 'spatialorder'
       SpatialOrder   = varargin{Iarg+1};
    case 'binsizespatial'
       BinSizeSpatial = varargin{Iarg+1};
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

[Ny,Nx] = size(ImageMatrix);

if (isempty(StartY)==1),
   % select StartY interactively
   disp('--- Open ds9 ---')
   % open FITS file in ds9 display   
   ds9_disp(ImageName,1);
   disp('--- Mark spatial position on image in which to start wavelength calibration');
   [Xspat,StartY] = ds9_getcoo(1,'image');
end

HalfBinSizeSpatial = 0.5.*(BinSizeSpatial - 1); 
VecStartY = [StartY-HalfBinSizeSpatial: 1 :StartY+HalfBinSizeSpatial];

% intensity as a function of dispersion position (uncalibrated
% wavelength=pixels) at the StartY spatial position.
IntensityDispDir_AtY = median(ImageMatrix(VecStartY,:),1).';
VecX                 = [1:1:Nx].';


% If preliminary line lists exist call: wavecalib_1d_auto.m
% If calibrated spectrm exist call: wavecalib_1d_xcorr.m
% else call: wavecalib_1d_int.m, when exit call: wavecalib_1d_auto.m

[SelWaveCalib,WaveCalibOut] = wavecalib_1d(IntensityDispDir_AtY,'Ignore','y',varargin);


% run the same wave calib in additional spatial positions (Y).


% Output is a matrix of wavelength per each pixel position.
MatW=interp2(Xw,Yw,W,MatX,MatY,InterpMethod,InterpMethod);
