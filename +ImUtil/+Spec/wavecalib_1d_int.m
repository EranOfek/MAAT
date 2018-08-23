function [SelWaveCalib,WaveCalibOut,RMS,ParWaveCalib]=wavecalib_1d_int(IntensityDispDir_AtY,varargin)
%-----------------------------------------------------------------------
% waveclib_1d_int function                                       ImSpec
% Description: Perform interactive wavelength calibration for 1-D spectrum.
% Input  : - Column vector containing 1-D spectrum [Intensity].
%            If two columns the [dispersion position(pix), Intensity].
%          *
%            'WaveCalib'       - Wavelength calibration matrix, containing
%                                [dispersion_position(X), Wavelength(A)],
%                                default is no calibration matrix (i.e., [];
%                                In this case the wavelength calibration
%                                will be preform manually).
%            'SigmaClip'       - Sigma clipping [-Number_sigma_low, +Number_sigma_high],
%                                for rejection of mis-identified lines,
%                                default is [-2, 2].
%            'PeakSearch'      - Semiwidth of line peak region size,
%                                default is 10 pixels.
%            'CenterMethod'    - Line peak centering method {'max'|'wmean'},
%                                default is 'stirling4'.
%            'WaveOrder'       - Polynomial order for fit in wavelength
%                                direction, default is [].
%                                If empty matrix (i.e., []), then prompt
%                                user for polynomial order.
%            'PlotLines'       - Plot identified arc lines {'y'|'n'},
%                                default is 'y'.
%            'PlotFit'         - Preform fit interactively, and plot best
%                                fit for wavelength calibration,
%                                {'y'|'n'}, default is 'y'.
%            'VerifyKey'       - Verify incorrect keywords {'y'|'n'},
%                                default is 'y'.
% Output : - List of the peak wavelength for the selected line peaks
%            [X, Wavelength, ErrX, Y, FlagMax].
%          - Calibration matrix for each pixel [X, Wavelength].
%          - RMS of final wavelength calibration fit.
%          - Parameters of wavelength calibration fit (from fit_lin.m).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Mar 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%-----------------------------------------------------------------------
if (nargin==0),
   error('Illegal number of input arguments');
end

Nx = size(IntensityDispDir_AtY,1);

if (size(IntensityDispDir_AtY,2)==1),
   VecX = [1:1:Nx].';
else
   % assume IntensityDispDir_AtY contains two columns
   VecX = IntensityDispDir_AtY(:,1);
   IntensityDispDir_AtY = IntensityDispDir_AtY(:,1);
end

% set default parameters
DefV.WaveCalib        = [];
DefV.SigmaClip        = [-2, 2];
DefV.PeakSearch       = 10;
DefV.CenterMethod     = 'stirling4';
DefV.WaveOrder        = 3;
DefV.PlotLines        = 'y';
DefV.PlotFit          = 'y';
DefV.VerifyKey        = 'y';

InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

% intensity as a function of dispersion position (uncalibrated
% wavelength=pixels) at the StartY spatial position.



%-------------------------------------
%--- Manual wavelength calibration ---
%-------------------------------------
figure;
plot(VecX,IntensityDispDir_AtY);
disp('--- Manual wavelength calibration ---');

if (isempty(InPar.WaveCalib)==1),
   % Input WaveCalib not exist
   SelWaveCalib = zeros(0,5);  % [X, UserWavelength, ErrX, Y, FlagMax]
else
   % Input WaveCalib exist
   Nw  = size(InPar.WaveCalib,1);
   Nc  = size(InPar.WaveCalib,2);
   SelWaveCalib = [InPar.WaveCalib, zeros(Nw,5-Nc)];
end
UWI    = size(SelWaveCalib,1);
Button = 1;
while Button==1,
   if (UWI>=2),   
      % on the fly linear fit for wavelength calibration
      ParWaveCalib=polyfit(SelWaveCalib(:,1),SelWaveCalib(:,2),1);
      Steps = 1000;   % A  - presenting solution at steps...
      Xs = [floor(min(VecX)./Steps).*Steps: Steps : ceil(max(VecX)./Steps).*Steps].';
      Ys = polyval(ParWaveCalib,Xs);
      disp([Xs Ys])
   end

   zoom on;
   Ans = input('Zoom in/out - press any key to continue','s')
   disp(sprintf('left click for selection, right click to exit interactive wavelength calibration'))
   [X0,Y,Button] = ginput(1);

   if (Button==1),
      % do not exit
      % find peak near mouse click
      [MaxX,MaxY,ErrX,FlagMax] = find_peak(VecX,IntensityDispDir_AtY,X0,...
                                           InPar.CenterMethod,InPar.PeakSearch);

      if (abs(FlagMax)==0),
         % maximum is OK
         UserWave = input('Enter corresponding wavelength (NaN to reject): ','s');
         UserWave = sscanf(UserWave,'%f');
         if (isnan(UserWave)==0),
            % Save maximum
 	    UWI = UWI + 1;
            % [X, UserWavelength, ErrX, Y, FlagMax]
            SelWaveCalib(UWI,:) = [MaxX, UserWave, ErrX, MaxY, FlagMax];
            % print to screen
            disp(sprintf('  X=%7.2f+/-%7.2f   W=%7.2f    Flag=%2d',MaxX,ErrX,UserWave,FlagMax))
         end
      else
         % maximum may be wrong
         disp('* cannot find maximum in search region');
      end
   end
end


%-----------------------------------------------------------
%--- construct the polynomial for wavelength calibration ---
%-----------------------------------------------------------
Nwc = size(SelWaveCalib,1);
OrigWaveOrder = InPar.WaveOrder;
InPar.WaveOrder = min(Nwc-1,InPar.WaveOrder);  % resetting WaveOrder
disp(sprintf('Reset WaveOrder from %d to %d',OrigWaveOrder,InPar.WaveOrder))



%----------------------------------------
%--- final wavelength calibration fit ---
%----------------------------------------
switch lower(InPar.PlotFit)
 case 'y'
    % interactive
    MaxIter = 1;
    PlotOption = 1;
    [XY2,XY1,XY0,ParWaveCalib,Resid,IU,INU,GraphicH] = plot_rm(SelWaveCalib(:,1),SelWaveCalib(:,2),'b.','polyfit_sc',5,InPar.WaveOrder,InPar.SigmaClip,MaxIter,PlotOption);

    %Res=plot_int({SelWaveCalib(:,1),SelWaveCalib(:,2)   
 %   ,IndRM,CallFun,FunPar,FunParInd,FunBehav)
    
 case 'n'
    % non interactive
    [ParWaveCalib,Resid]      = polyfit_sc(SelWaveCalib(:,1),SelWaveCalib(:,2),InPar.WaveOrder,InPar.SigmaClip);

 otherwise
    error('Unknown PlotFit option');
end

RMS=std(Resid);

%--- Plot the spectrum and all identified lines ---
switch lower(InPar.PlotLines)
 case 'y'
    figure;
    plot(VecX,IntensityDispDir_AtY);
    hold on;
    plot(SelWaveCalib(:,1),   SelWaveCalib(:,4),   'ro');
    if (isempty(InPar.WaveCalib)==0),
       plot(   InPar.WaveCalib(:,1),   InPar.WaveCalib(:,4),      'rs','MarkerSize',12);
    end
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



