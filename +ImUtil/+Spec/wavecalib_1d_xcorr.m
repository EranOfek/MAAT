function [SelWaveCalib,WaveCalibOut,RMS,ParWaveCalib]=wavecalib_1d_xcorr(IntensityDispDir_AtY,ArcTemplate,varargin)
%-----------------------------------------------------------------------
% waveclib_1d_xcorr function                                     ImSpec
% Description: Perform automatic wavelength calibration for 1-D spectrum,
%              by cross corelating it with a template spectrum.
% Input  : - Vector containing 1-D spectrum [Intensity].
%            If two columns the [dispersion position(pix), Intensity].
%          - Arc template to cross correlate with spectrum,
%            [wavelength(Ang), Intensity].
%            If string then, a mat file name containing the arc template
%            matrix.



%          *

%            'WaveCalib'       - Wavelength calibration matrix, containing
%                                [dispersion_position(X), Wavelength(A)],
%                                default is no calibration matrix (i.e., [];
%                                In this case the wavelength calibration
%                                will be preform manually).
%            'SigmaClip'       - Sigma clipping [-Number_sigma_low, +Number_sigma_high],
%                                for rejection of mis-identified lines,
%                                default is [-2, 2].
%            'Scale'           - Optional scale in Ang/pix of the input
%                                1-D spectrum to cross calibrate.
%                                By default the scale is set to the same
%                                scale as the ArcTemplate.
%            'RefWave'         - Approximate reference pixel and wavelength
%                                in the input spectrum,
%                                [Ref_Pixel(pixel), Ref_wavelength(Ang)].
%                                By default taking: [1, ArcTemplate(1,1)].
%            'XCorrWindow'     - The half window size [Ang] in which
%                                to cross-correlate the spectrum and template.
%                                Default is 100.
%            'FindShift'       - Method for finding shift between template
%                                and spectrum {'CCF','Interactive'},
%                                default is 'CCF'.
%            'PeakSearch'      - Semiwidth of line peak region size,
%                                default is 10 pixels.
%            'CenterMethod'    - Line peak centering method {'max'|'wmean'},
%                                default is 'wmean'.
%            'WaveOrder'       - Polynomial order for fit in wavelength
%                                direction, default is 3.
%            'PlotCCF'         - Plot CCF {'y'|'n'}, default is 'y'.
%            'PlotLines'       - Plot identified arc lines {'y'|'n'},
%                                default is 'y'.
%            'PlotFit'         - Plot best fit for wavelength calibration,
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
import Util.fit.*

if (nargin==0),
   error('Illegal number of input arguments');
end

if (ischar(ArcTemplate)==1),
   eval(sprintf('load %s',ArcTemplate));
   eval(sprintf('ArcTemplate=%s;',ArcTemplate));
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
WaveCalib        = [];
SigmaClip        = [-2, 2];
Scale            = ScaleTemplate;
RefWave          = [1, ArcTemplate(1,1)];
XCorrWindow      = 100;
FindShift        = 'CCF';
PeakSearch       = 10;
CenterMethod     = 'wmean';
WaveOrder        = 3;
PlotCCF          = 'y';
PlotLines        = 'y';
PlotFit          = 'y';
Ignore           = 'n';

InPar = set_varargin_keyval(DefV,'y','use',varargin{:});


% Interpolate ArcTemplate such that it will be evenly spaced with
% the smallest step.
MinStep = min(diff(ArcTemplate(:,1)));
MaxStep = max(diff(ArcTemplate(:,1)));
if (MinStep==MaxStep),
   % do nothing: ArcTemplate is already evenly spaced
else
   MinWave = min(ArcTemplate(:,1));
   MaxWave = max(ArcTemplate(:,1));
   Wave    = [MinWave:MinStep:MaxWave]';
   ArcTemplate = [Wave, interp1(ArcTemplate(:,1),ArcTemplate(:,2),Wave,InPar.ArcInterp)];
end

got here
....................................................................




ScaleTemplate    = abs(ArcTemplate(2,1) - ArcTemplate(1,1));


% intensity as a function of dispersion position (uncalibrated
% wavelength=pixels) at the StartY spatial position.



%-----------------------------------------------------------
%--- construct the polynomial for wavelength calibration ---
%-----------------------------------------------------------
WaveFun{1} = '1';
MapFun{1}  = [1];
for Io=2:1:WaveOrder+1,
   WaveFun{Io} = sprintf('x.^%d',Io-1);
   MapFun{Io}  = [1];
end


%----------------------------------------------
%--- Automatic xcorr wavelength calibration ---
%----------------------------------------------
disp('--- XCorr wavelength calibration ---');




switch lower(FindShift)
 case 'ccf'
    % find shift between template and spectrum by cross-correlation
    N_AT = size(ArcTemplate,1);
    % ArcTemplate  [Wave, Intensity]
    Spec = [(RefWave(1)-VecX).*ScaleTemplate./Scale + RefWave(2), IntensityDispDir_AtY]
????
    CCF  = ccf([1:1:N_AT].'.*Scale./ScaleTemplate, ArcTemplate(:,2)],...
               [VecX, IntensityDispDir_AtY], ...
	        1, 'normal', 'n', XCorrWindow);

    swith PlotCCF
     case 'y'
        figure;
        plot(CCF(:,1),CCF(:,2));
     case 'n'
        % do nothing
     otherwise
        error('Unknown PlotCCF option');
    end

    [MaxCorr,MaxInd] = max(CCF(:,2));
    WaveShift        = CCF(MaxInd,1).*Scale;   % [Ang]



...

 case {'int','interactive'}
    error('Not implemented yet');

 otherwise
    error('Unknown FindShift option');
end




%----------------------------------------
%--- final wavelength calibration fit ---
%----------------------------------------
switch PlotFit
 case 'y'
    PlotOptions = [0 1 0];
 case 'n'
    PlotOptions = [0 0 0];
 otherwise
    error('Unknown PlotFit option');
end

% First iteration fit - no sigma cliping
[ParWaveCalib,PE,Chi2,Dof,Resid,Cov,IndSigClip] = fit_lin(WaveFun,SelWaveCalib(:,2),0.1,SelWaveCalib(:,1),MapFun,'SVD','Plot',PlotOptions);
title('first iteration wavelength calibration fit');

% Second iteration fit - apply sigma cliping
Isc  = find(Resid>=-abs(SigmaClip(1)).*std(Resid) & Resid<=abs(SigmaClip(2)).*std(Resid));
Insc = find(Resid<-abs(SigmaClip(1)).*std(Resid) | Resid>abs(SigmaClip(2)).*std(Resid));

[ParWaveCalib,PE,Chi2,Dof,Resid,Cov,IndSigClip] = fit_lin(WaveFun,SelWaveCalib(Isc,2),0.1,SelWaveCalib(Isc,1),MapFun,'SVD','Plot',PlotOptions);
title('second iteration wavelength calibration fit');
RMS=std(Resid);

%--- Plot the spectrum and all identified lines ---
switch PlotLines
 case 'y'
    figure;
    plot(VecX,IntensityDispDir_AtY);
    hold on;
    plot(SelWaveCalib(:,1),   SelWaveCalib(:,4),   'ro');
    plot(   WaveCalib(:,1),   WaveCalib(:,4),      'rs','MarkerSize',12);
    plot(SelWaveCalib(Insc,1),SelWaveCalib(Insc,4),'cx','MarkerSize',10);
    hold off;
 case 'n'
    % do nothing
 otherwise
    error('Unknown PlotLines option');
end


CalibW = lin_fun(WaveFun,ParWaveCalib,VecX,MapFun);
WaveCalibOut = [VecX, CalibW];
