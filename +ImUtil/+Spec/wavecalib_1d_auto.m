function [SelWaveCalib,WaveCalibOut,RMS,ParWaveCalib]=wavecalib_1d_auto(IntensityDispDir_AtY,varargin)
%-----------------------------------------------------------------------
% waveclib_1d_auto function                                      ImSpec
% Description: Perform automatic wavelength calibration for 1-D spectrum.
% Input  : - Vector containing 1-D spectrum [Intensity].
%            If two columns the [dispersion position(pix), Intensity].
%          *
%            'WaveCalib'       - Wavelength calibration matrix, containing
%                                [dispersion_position(X), Wavelength(A)],
%                                default is no calibration matrix (i.e., [];
%                                In this case the wavelength calibration
%                                will be preform manually).
%            'ArcType'         - Arc lines type (see get_arclines.m
%                                for more details).
%                                Options are: {'Ar'|'Cd'|'Hg'|'Ne'|'Zn'|
%                                              'SA'|'SS'}, default is 'SS'.
%            'SigmaClip'       - Sigma clipping [-Number_sigma_low, +Number_sigma_high],
%                                for rejection of mis-identified lines,
%                                default is [-2, 2].
%            'PeakSearch'      - Semiwidth of line peak region size,
%                                default is 10 pixels.
%            'CenterMethod'    - Line peak centering method {'max'|'wmean'|'stirling4'},
%                                default is 'wmean'.
%            'WaveOrder'       - Polynomial order for fit in wavelength
%                                direction, default is 3.
%                                If empty matrix (i.e., []), then prompt
%                                user for polynomial order.
%            'PlotLines'       - Plot identified arc lines {'y'|'n'},
%                                default is 'y'.
%            'PlotFit'         - Preform fit interactively, and plot best
%                                fit for wavelength calibration,
%                                {'y'|'n'}, default is 'y'.
%            'Ignore'          - Ignore incorrect keywords {'y'|'n'},
%                                default is 'n'.
%                                If 'y', then should be first keyword.
% Output : - List of the peak wavelength for the selected line peaks
%            [X, Wavelength, ErrX, Y, FlagMax].
%          - Calibration matrix for each pixel [X, Wavelength].
%          - RMS of final wavelength calibration fit.
%          - Parameters of wavelength calibration fit (from fit_lin.m).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                     March 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------
Thresh = 0.1;

if (nargin==0)
   error('Illegal number of input arguments');
end

Nx = size(IntensityDispDir_AtY,1);

if (size(IntensityDispDir_AtY,2)==1)
   VecX = [1:1:Nx].';
else
   % assume IntensityDispDir_AtY contains two columns
   VecX = IntensityDispDir_AtY(:,1);
   IntensityDispDir_AtY = IntensityDispDir_AtY(:,1);
end

% set default parameters
WaveCalib        = [];
ArcType          = 'SS';
SigmaClip        = [-2, 2];
PeakSearch       = 10;
CenterMethod     = 'stirling4';
WaveOrder        = 3;
PlotLines        = 'y';
PlotFit          = 'y';
Ignore           = 'n';

Narg = length(varargin);
for Iarg=1:2:Narg
   switch lower(varargin{Iarg})
    case 'wavecalib'
       WaveCalib      = varargin{Iarg+1};
    case 'arctype'
       ArcType        = varargin{Iarg+1};
    case 'sigmaclip'
       SigmaClip      = varargin{Iarg+1};
    case 'peaksearch'
       SPeakSearch    = varargin{Iarg+1};
    case 'centermethod'
       CenterMethod   = varargin{Iarg+1};
    case 'waveorder'
       WaveOrder      = varargin{Iarg+1};
    case 'plotlines'
       PlotLines      = varargin{Iarg+1};
    case 'plotfit'
       PlotFit        = varargin{Iarg+1};
    case 'ignore'
       Ignore         = varargin{Iarg+1};
    otherwise
       switch lower(Ignore)
        case 'n'
           error('Unknown keyword option');
        case 'y'
           % do nothing
        otherwise
	   error('Unkown Ignore option');
       end
   end
end


% intensity as a function of dispersion position (uncalibrated
% wavelength=pixels) at the StartY spatial position.



% get arclines
ArcLines = get_arclines(ArcType);

%----------------------------------------
%--- Automatic wavelength calibration ---
%----------------------------------------
disp('--- Automatic wavelength calibration ---');

if (isempty(WaveCalib)==1)
   % Input WaveCalib not exist
   SelWaveCalib = zeros(0,5);  % [X, UserWavelength, ErrX, Y, FlagMax]
   error('wavecalib_1d_auto cannot run without initial guess');
else
   % Input WaveCalib exist
   Nw           = size(WaveCalib,1);
   Nc           = size(WaveCalib,2);
   SelWaveCalib = [WaveCalib, zeros(Nw,5-Nc)];

   Narc = size(ArcLines,1);
   MI = xyxymatch([zeros(Nw,1),   SelWaveCalib(:,2)],...
                  [zeros(Narc,1), ArcLines],...
                  [1 2],[1 2],Thresh,Thresh,'plane');
   In = find(isnan(MI(:,2))==0);
   FoundI  = MI(In,2);   % indices of pre-identified arc lines
   ArcFlag = zeros(Narc,1);
   ArcFlag(FoundI) = 1;  

end
UWI    = size(SelWaveCalib,1);
MouseB = 1;

AllFound = 0;
while AllFound==0
   if (UWI>=2)
      % on the fly linear fit for wavelength calibration
      ParWaveCalib=polyfit(SelWaveCalib(:,1),SelWaveCalib(:,2),1);
      Steps = 1000;   % A  - presenting solution at steps...
      Xs = [floor(min(VecX)./Steps).*Steps: Steps : ceil(max(VecX)./Steps).*Steps].';
      Ys = polyval(ParWaveCalib,Xs);
      %disp([Xs Ys])
   end

   [NearestW,IndW,DeltaW]=nearest_unflag(SelWaveCalib(1,2),ArcLines,ArcFlag);
   ArcFlag(IndW) = 1;    % flag arc line as identified

   if (isempty(NearestW)==1)
      AllFound = 1;   % All lines identified - exit
   else
      % do not exit
      % find nearest line peak

      % convert NearestW to position along dispersion axis
      Roots = roots([ParWaveCalib(1:end-1),ParWaveCalib(end)-NearestW]);
      Iroot = find(Roots>min(VecX) & Roots<max(VecX) & isreal(Roots)==1);
      if (length(Iroot)~=1),
         % skip - less or more than one solution
      else
         X0 = Roots(Iroot);
         [MaxX,MaxY,ErrX,FlagMax] = find_peak(VecX,IntensityDispDir_AtY,X0,CenterMethod,PeakSearch,'highest',PeakSearch);

         if (abs(FlagMax)==0),
            % maximum is OK

            % Save maximum
            UWI = UWI + 1;
            % [X, UserWavelength, ErrX, Y, FlagMax]
            SelWaveCalib(UWI,:) = [MaxX, NearestW, ErrX, MaxY, FlagMax];
            % print to screen
            disp(sprintf('  X=%7.2f+/-%7.2f   W=%7.2f    Flag=%2d',MaxX,ErrX,NearestW,FlagMax))
         end
      end
   end
end


%----------------------------------------
%--- final wavelength calibration fit ---
%----------------------------------------
switch PlotFit
 case 'y'
    % interactive
    MaxIter = 1;
    PlotOption = 1;
    [XY2,XY1,XY0,ParWaveCalib,Resid,IU,INU,GraphicH] = plot_rm(SelWaveCalib(:,1),SelWaveCalib(:,2),'b.','polyfit_sc',5,WaveOrder,SigmaClip,MaxIter,PlotOption);

 case 'n'
    % non interactive
    [ParWaveCalib,Resid]      = polyfit_sc(SelWaveCalib(:,1),SelWaveCalib(:,2),WaveOrder,SigmaClip);

 otherwise
    error('Unknown PlotFit option');
end

RMS=std(Resid);

%--- Plot the spectrum and all identified lines ---
switch PlotLines
 case 'y'
    figure;
    plot(VecX,IntensityDispDir_AtY);
    hold on;
    plot(SelWaveCalib(:,1),   SelWaveCalib(:,4),   'ro');
    plot(   WaveCalib(:,1),   WaveCalib(:,4),      'rs','MarkerSize',12);
    plot(SelWaveCalib(:,1),SelWaveCalib(:,4),'cx','MarkerSize',10);
    hold off;
 case 'n'
    % do nothing
 otherwise
    error('Unknown PlotLines option');
end


%CalibW = lin_fun(WaveFun,ParWaveCalib,VecX,MapFun);
CalibW = polyval(ParWaveCalib,VecX);
WaveCalibOut = [VecX, CalibW];
