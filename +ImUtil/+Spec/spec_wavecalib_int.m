function [WaveCalib,SelectedLines,Poly]=spec_wavecalib_int(IntensityDispDir_AtY,varargin)
%------------------------------------------------------------------------------
% spec_wavecalib_int function                                           ImSpec
% Description: Perform interactive wavelength calibration for 1-D spectrum.
% Input  : - Column vector containing 1-D spectrum [Intensity].
%            If two columns the [dispersion position(pix), Intensity].
%          * Pairs of ...,key,val,... input arguments.
%            The following keywords are available:
%            'PeakSearch'   - Semi width of the region in which to search
%                             for the line peak around the user clicked
%                             position. Default is 10 pixels.
%            'CenterMethod' - Line peak centering method.
%                             See find_peak.m for options.
%                             Default is 'stirling4'.
%            'Deg'          - Degree of default polynomial fit.
%                             Default is 2.
%            'FunPar'       - Cell array of additional parameters to pass
%                             to fitgenpoly.m.
%                             Default is {'Method','StdP','Clip',[2 2],
%                             'MaxIter',1,'Mean','Median','Plot','none'};
%            'IntFFit'      - Perform the final fit interactively {'y'|'n'}.
%                             Default is 'y'.
%            'VerifyKey'    - Verify incorrect keywords {'y'|'n'},
%                             default is 'y'.
% Output : - A structure containing the wavelength calibration for each pixel.
%            The following fields are available:
%            .X       - Vector of X (pix) positions.
%            .Int     - Vector of the input spectrum intensity.
%            .Wave    - Vector of wavelength [Ang] corresponding to each pixel
%                       in .X.
%            .WaveErr - Error in wavelngth [Ang] corresponding to each pixel
%                       in .X.
%            .RMS     - RMS of the best fit solution.
%          - Structure array of the lines selected by the user.
%            Each element contains the following fields are available:
%            .ClickX   - The X position clicked by the user.
%            .ClickY   - The Y position clicked by the user.
%            .MaxX     - The best X position of the line peak.
%            .MaxY     - The best Y position of the line peak.
%            .LineWave - The line wavelength entered by the user.
%          - The strucure of the parameters of the best fit polynomial.
%            See fitgenpoly.m for details.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Mar 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [WaveCalib,SelectedLines,Poly]=spec_wavecalib_int(IntY);
% S=fitsread('lred0178.fits');
% IntY=(median(S(200:210,:)))';

%------------------------------------------------------------------------------
NpMax = 3;
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
DefV.PeakSearch       = 10;
DefV.CenterMethod     = 'stirling4';
DefV.IntFFit          = 'y'; % perform final fit interactively
DefV.Deg              = 2;
DefV.FunPar           = {'Method','StdP','Clip',[2 2],'MaxIter',1,'Mean','Median'};
DefV.VerifyKey        = 'y';

InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

% intensity as a function of dispersion position (uncalibrated
% wavelength=pixels) at the StartY spatial position.



%-------------------------------------
%--- Manual wavelength calibration ---
%-------------------------------------
Hfig   = figure;
Hlines = plot(VecX,IntensityDispDir_AtY);
Hgca   = gca;
H      = xlabel('Wavelength [pix]');
set(H,'FontSize',16);
H      = ylabel('Intensity');
set(H,'FontSize',16);

MarkerLength = max(IntensityDispDir_AtY).*0.02;  % length of line marker

%if (isempty(InPar.WaveCalib)==1),
%   % Input WaveCalib not exist
%   SelWaveCalib = zeros(0,5);  % [X, UserWavelength, ErrX, Y, FlagMax]
%else
%   % Input WaveCalib exist
%   Nw  = size(InPar.WaveCalib,1);
%   Nc  = size(InPar.WaveCalib,2);
%   SelWaveCalib = [InPar.WaveCalib, zeros(Nw,5-Nc)];
%end
%UWI    = size(SelWaveCalib,1);

fprintf('Start interactive wavelength calibration\n');
fprintf('----------------------------------------\n');
plot_wavecalib_int_menu

PolyTry = [];
SelectedLines = [];
QuitInt = false;
Iline   = 0;
while (~QuitInt);
   fprintf('Click keyboard key to continue - (''m'' for menu)\n');
   IntRes=plot_int1(Hfig,'key');
   switch lower(IntRes.Key)
    case 'q'
       % Quit
       QuitInt = true;
    case {'m','?'}
       % show menu
       plot_wavecalib_int_menu;
    case {'l'}
       % show a list of selected lines
       Nl = length(SelectedLines);
       Np = length(PolyTry);
       Resid = zeros(Nl,NpMax).*NaN;
       for Ip=1:1:Np,
  	  Resid(:,Ip) = PolyTry(Ip).Resid;
       end

       fprintf('     X_peak  Wavelength        Residuals\n');
       fprintf('  ---------- ----------    P(1)    P(2)    P(3)\n');
       for Inl=1:1:Nl,
           % note that this will work properly only with NpMax=3
 	   fprintf('  %10.3f %10.3f   %7.3f %7.3f %7.3f\n',...
			 SelectedLines(Inl).MaxX,...
			 SelectedLines(Inl).LineWave,...
                         Resid(Inl,:));
       end


    case 'i'
       % attempt to identify lines automatically.
       fprintf('i NOT implemented');
       WaveCalibLinesPar = {'ContrastRMS',8,'MaxLineShift',10,...
                            'FitPoly','n','PlotMS','n','PlotFit',n'};
       % 
       
%PolyP = Poly.FunOut{1};

%WaveCalib.X       = VecX;
%WaveCalib.Wave    = polyval(PolyP.Par,VecX);
%%[Y,DelY]          = polyconf_cov(PolyP.Par,VecX,PolyP.Cov);
%%WaveCalib.WaveErr = DelY;
%%WaveCalib.RMS     = PolyP.RMS;

%SpecIn = [WaveCalib.Wave, IntensityDispDir_AtY];
%InPar.Template = 'Ar';

%       [PolyML,MatchedListML] = spec_wavecalib_lines(SpecIn,InPar.Template,WaveCalibLinesPar{:});
%MatchedListML

    case 's'
       % select line follows by prompting for line wavelength
       fprintf('Use left mouse click to select line center\n');
       IntRes   = plot_int1(Hfig,'mouse');
       WaveStr  = input(' Enter spectral line wavelength [A] : ','s');
       LineWave = str2num(WaveStr); 

       [MaxX,MaxY,ErrX,FlagMax] = find_peak(VecX,IntensityDispDir_AtY,...
                                            IntRes.Pos(1),...
                                            InPar.CenterMethod,...
                                            InPar.PeakSearch);
       if (abs(FlagMax)==0),
          % max is ok
          % mark line center on plot
          hold on;
          Hsel  = plot([MaxX MaxX],MaxY+[0 MarkerLength],'r-');
          Iline = Iline + 1;
          SelectedLines(Iline).ClickX   = IntRes.Pos(1);
          SelectedLines(Iline).ClickY   = IntRes.Pos(2);
          SelectedLines(Iline).MaxX     = MaxX;
          SelectedLines(Iline).MaxY     = MaxY;
          SelectedLines(Iline).LineWave = LineWave;
          SelectedLines(Iline).Hsel     = Hsel;

          % on the fly polynomial fit
          PolyTry = fit_various_wavecalib([SelectedLines.MaxX].',...
                                          [SelectedLines.LineWave].',NpMax);

       else
          fprintf('Can not find line peak in range\n');
       end
    case 'd'
       fprintf('Use left mouse click to delete selected line\n');
       IntRes   = plot_int1(Hfig,'mouse');
       [MinDist,MinInd] = min(abs(IntRes.Pos(1)-[SelectedLines(:).MaxX]));
       % delete line marker handle from plot
       delete(SelectedLines(MinInd).Hsel);
       Flag          = boolean(zeros(Iline,1));
       Flag(MinInd)  = 0;
       SelectedLines = SelectedLines(Flag);
       Iline         = length(SelectedLines);

       % on the fly polynomial fit
       PolyTry = fit_various_wavecalib([SelectedLines.MaxX].',...
                                       [SelectedLines.LineWave].',NpMax);

    otherwise
       % show menu
       plot_wavecalib_int_menu;
   end
end      

% close spectrum figure
delete(Hfig);

% Perform final wavecalib fitting
switch lower(InPar.IntFFit)
 case 'y'
    % interactive final fitting
    figure;
    Hfig = gcf;
    Hf = plot([SelectedLines.MaxX].',...
              [SelectedLines.LineWave].','ko');
    H  = xlabel('Wavelength [pix]');
    set(H,'FontSize',16);
    H  = ylabel('Wavelength [Ang]');
    set(H,'FontSize',16);

    Poly = plot_int(Hfig,[],'fitgenpoly',...
                    {1, InPar.Deg, InPar.FunPar{:},'Plot','fitonly'},...
                    'FunParInd',2,'FunBehav','i','DispFit','y');
    waitfor(Hfig,'KeyPressFcn','');
 case 'n'
    Poly = Util.fit.fitgenpoly([SelectedLines.MaxX].',...
                      [SelectedLines.LineWave].',...
                      1,InPar.Deg,InPar.FunPar{:});
 otherwise
    error('Unknown IntFFit option');
end


PolyP = Poly.FunOut{1};

% set the Wavecalib matrix
WaveCalib.X       = VecX;
WaveCalib.Int     = IntensityDispDir_AtY;  % input spectrum instensity
WaveCalib.Wave    = polyval(PolyP.Par,VecX);
[Y,DelY]          = polyconf_cov(PolyP.Par,VecX,PolyP.Cov);
WaveCalib.WaveErr = DelY;
WaveCalib.RMS     = PolyP.RMS;


return


%----------------------------------------------------------------------
function Poly=fit_various_wavecalib(X,Y,NpMax,varargin)
%----------------------------------------------------------------------
%NpMax  = 3;
Np     = min(NpMax,length(X)-1);
Err    = 1;

Poly = [];
if (Np>0),
   fprintf('On the fly polynomial fits using %d points:\n',length(X));
   clear Poly;
   for Ip=1:1:Np,
      Poly(Ip) = Util.fit.fitgenpoly(X,Y,Err,Ip,varargin{:});
      fprintf('   Polyfit deg=%d, Slope=%f, RMS=%f\n',...
              length(Poly(Ip).Par)-1,...
	      Poly(Ip).Par(end-1),...
	      Poly(Ip).RMS);
   end
end

return


%----------------------------------------------------------------------
function plot_wavecalib_int_menu
%----------------------------------------------------------------------
fprintf('Menu:\n');
fprintf('   Click ''s'' for selection, follows by left click on line position and entering line wavelength\n');
fprintf('   Click ''d'' to delete selection, follows by left click on line position\n');
fprintf('   Click ''i'' to identify lines automatically based on initial solution\n');
fprintf('   Click ''l'' to show a list of selected lines\n');
fprintf('   Click ''q'' to quit\n');
fprintf('   Click ''m'' to show this menu\n');
fprintf('\n');

return

