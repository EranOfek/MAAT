function [PolyFit,AllLines]=spec_wavecalib_xcorr(Spec,ArcTemplate,varargin)
%-----------------------------------------------------------------------
% waveclib_1d_xcorr function                                     ImSpec
% Description: Perform automatic wavelength calibration for 1-D spectrum,
%              by cross corelating it with a template spectrum.
%              The cross-correlation is done using various trial scales,
%              so the program optimize for both the scale and shift.
%              After the cross-correlation, refine the wavelength
%              calibration by line matching.
% Input  : - A column vector containing 1-D spectrum [Intensity].
%            If two columns the [dispersion position(pix), Intensity].
%            Note that the dispersion position must be equaly spaced.
%          - Arc template to cross correlate with spectrum,
%            [wavelength(Ang), Intensity].
%            If string then one of the following arcs available via
%            the spec_get_arc.m function:
%            {'SkyLow'|'SkyHigh'|'Ar'|'Cd'|'Hg'|'Ne'|'Zn'}.
%            For a full list and details see spec_get_arc.m
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Back'   - Method by which to subtract background from
%                       the spectrum prior to the cross correlation.
%                       See subtract_back1d.m for options.
%                       Default is 'median'.
%            'BackPar' - Additional background subtraction parameters
%                       to pass to subtract_back1d.m. Default is NaN.
%            'ScaleVec' - Vector of scales to test.
%                       Default is logspace(-1,1,3000).'.
%            'RefineXC' - Refine the cross correlation using a finder grid
%                       around best solution {'y'|'n'}. Default is 'y'.
%            'Refinment' - The refinment scale step size measured in
%                       absolute units.
%                       Default is 0.0001.
%            'InterpMethod' - Interpolation method. Default is 'linear'.
%            'MaxLineShift' - Maximum allowed shift between predicted and
%                       calculated line positions in the first iteration
%                       of line matching. Default is 10 A.
%            'MaxLineShiftI' - Maximum allowed shift between predicted and
%                       calculated line positions in the next iterations
%                       of line matching. Default is 3 A.
%            'PlotXC'   - Plot cross correlation results {true|false}.
%                       Default is false.
%            'FitLinesNiter' - Number of line matching iterations.
%                       Default is 2.
%            'PlotMS' - Plot spectra with matched lines {'n'|'y'}.
%                       Default is 'n'.
%            'PlotFit' - Plot best fit line matching {'n'|'y'}.
%                       Default is 'n'.
%            'Deg'     - Polynomial degree for the first iteartion of
%                       the line matching fit. Default is 5.
%            'DegI'    - Polynomial degree to the residuals fitting of
%                       the bext iterations. Default is 1.
%            'Interactive' - Interactive line matching mode {'n'|'y'}.
%                       Default is 'n'.
%            'PolyFitPar' - Additional parameters for the polynomial fit.
%                       Default is
%                       {'Algo','chol','NormX','y','MaxIter',1,'Method','StdP','Mean','median','Clip',[2 2]}.
% Output : - A structure containing the best fit parameters with the
%            following fields:
%            .ArcScale  - The scale of the original Arc template [A/pix].
%            .BestShift - Best shift.
%            .BestScale1- Best scale in units of the Arc scale.
%            .BestScale - Best scale of spectrum  [A/pix].
%            .BestPix1  - Wavelength [A] of first pixel in spectrum.
%            .xcSpecWave- Vector of wavelength for each pixel in spectrum
%                         based on the initial cross-correlation.
%            .SpecWave  - Vector of wavelength for each pixel in spectrum
%                         based on the final line matching.
%            .Poly      - Best fit polynomila parameters of the line
%                         matching. See fields description in fitgenpoly.m.
%          - Structure array of all matched lines in spectrum and
%            template. The structure contains the following fields:
%            .PossMatch - Number of possible matches in search window.
%            .TempWave  - Line wavelength in template [A].
%            .SpecWave  - Line wavelength in spectrum [A].
%            .IndSM     - Index of line in template lines list.
%          - ArcTemplate spectrum
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Fit,MatchedList]=spec_wavecalib_xcorr(dbsp_red,'Ne+Ar')
%          [Fit,MatchedList]=spec_wavecalib_xcorr(dbsp_red,'Ne+Ar','PlotFit','y')
% Reliable: 2
%-----------------------------------------------------------------------
if (nargin==0),
   error('Illegal number of input arguments');
end

if (~ischar(ArcTemplate)),
    % ArcTemplate is defined [Wave, Intensity]
else
   SpecArc = spec_get_arc(ArcTemplate);
   if (isempty(SpecArc)),
       error('Unknown ArcTemplate option');
   end
   ArcTemplate  = SpecArc.Spec;
   ListArcLines = SpecArc.Lines;
end

[Ni,Nj] = size(Spec);
if (Nj==1),
   SpecPix = (1:1:Ni).';
   SpecInt = Spec;
elseif (Nj>1),
   SpecPix = Spec(:,1);
   SpecInt = Spec(:,2);
else
   error('Spec input should be a column vector or two columns matrix');
end


%--- set default parameters ---
DefV.Back          = 'median';
DefV.BackPar       = NaN;
DefV.ScaleVec      = logspace(-1,1,3000).';
DefV.Nrange        = 4;
DefV.ScaleVec2     = (0.9:0.0001:1.1).';
DefV.MinCorr       = 0.7;
%DefV.Refinment     = 0.00001;
DefV.InterpMethod  = 'linear';
%DefV.ContrastRMS   = 8;
DefV.MaxLineShift  = 3;
%DefV.MaxLineShiftI = 3;   % maximum line shift after first iteration
DefV.PlotXC        = false;   % plot last XC 
DefV.RefineXC      = 'y';
DefV.FitLines      = 'y';
DefV.FitLinesNiter = 2;
% parameters for spec_wavecalib_lines.m
DefV.PlotMS       = 'n';   % plot spectra with matched lines
DefV.PlotFit      = 'n';   
DefV.Deg          = 5;
DefV.DegI         = 1;  % poly deg after first iteration
DefV.Interactive  = 'n';
DefV.PolyFitPar   = {'Algo','chol','NormX','y','MaxIter',1,'Method','StdP','Mean','median','Clip',[2 2]};

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});


%-----------------------------------------------------
%--- Subtract background from spectra before XCORR ---
%-----------------------------------------------------
SpecInt     = timeseries.subtract_back1d(SpecInt,InPar.Back,InPar.BackPar);
ArcTemplate = timeseries.subtract_back1d(ArcTemplate,InPar.Back,InPar.BackPar);

%Ns          = size(Spec,1);
%Nt          = size(ArcTemplate,1);

%------------------------------------
%--- XCORR as a function of scale ---
%------------------------------------
Info = xcorr_scale_shift([SpecPix,SpecInt],ArcTemplate,'ScaleVec',InPar.ScaleVec,'Back','medfilt','BackPar',500);
if (Info.BestCorr<InPar.MinCorr),
    % did not find good correlation
    PolyFit  = [];
    AllLines = [];
else

    Fit.xcSpecWave = Info.DataTemp;
    MaxCorrInd = Info.BestCol;
    SpecSec = [Fit.xcSpecWave, SpecInt];


    %--------------------------------------------------------
    %--- Run XCorr per section and find scale derivatives ---
    %--------------------------------------------------------
    Range = range(Fit.xcSpecWave)./InPar.Nrange;

    Min=min(Fit.xcSpecWave);
    max(Fit.xcSpecWave);

    AllLines = zeros(0,4);
    for Irange=1:1:InPar.Nrange,

       FlagRange       = Fit.xcSpecWave>(Min+(Irange-1).*Range) & Fit.xcSpecWave<(Min+Irange.*Range);
       FlagTemplate    = ArcTemplate(:,1)>(Min+(Irange-1).*Range) & ArcTemplate(:,1)<(Min+Irange.*Range);
       ScaleVec2 = Info.BestScale.*InPar.ScaleVec2; %(0.9:0.0001:1.1).';
       InfoAll(Irange) = xcorr_scale_shift([SpecPix(FlagRange),SpecInt(FlagRange)],ArcTemplate(FlagTemplate,:),'ScaleVec',ScaleVec2);
       SpecPixSection  = SpecPix(FlagRange);

    %    clf
    %    plot(InfoAll(Irange).DataTemp,SpecInt(FlagRange))
    %    hold on
    %    graph(ArcTemplate(FlagTemplate,:),'r-')
    %    input('any key','s');

        %--------------------------------------------------
        %--- Refine XCORR solution using specific lines ---
        %--------------------------------------------------

       [~,TmpLines]=spec_wavecalib_lines(SpecSec(FlagRange),ListArcLines,...
                                                'PlotMS',InPar.PlotMS,...
                                                'PlotFit',InPar.PlotFit,...
                                                'FitPoly','n',...
                                                'Interactive','n',...
                                                'MaxLineShift',InPar.MaxLineShift);

       PixPos   = interp1(InfoAll(Irange).DataTemp, SpecPixSection, [TmpLines.SpecWave].',InPar.InterpMethod);
       AllLines = [AllLines; [PixPos, [TmpLines.TempWave].', [TmpLines.SpecWave].', [TmpLines.PossMatch].']];

    end

    AllLines = AllLines(~isnan(AllLines(:,1)) & ~isnan(AllLines(:,2)) & AllLines(:,4)==1,:);
    polyfit(AllLines(:,1),AllLines(:,2),1)

    PolyFit = Util.fit.fitgenpoly(AllLines(:,1),AllLines(:,2),1,InPar.Deg,InPar.PolyFitPar{:});
    PolyFit.SpecWave = polyval(PolyFit.Par,SpecPix);

    %plot(AllLines(:,1),AllLines(:,2),'.')
    PolyFit.ArcTemplate = ArcTemplate;
    PolyFit.SpecPix     = SpecPix;
    PolyFit.SpecInt     = SpecInt;
end

