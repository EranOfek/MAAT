function [Poly,MatchedList]=spec_wavecalib_lines(Spec,Template,varargin)
%------------------------------------------------------------------------------
% spec_waveclib_lines function                                     ImSpec
% Description: Given a spectrum with an approximate wavelength solution,
%              identify lines in the spectrum and attempt to match them
%              with list of spectral lines in template spectra
%              (e.g., Arcs or sky lines). Then, refine the wavelength
%              calibration by fitting a function (e.g., polynomial) to
%              the lines in the spectrum to match the lines in the template.
% Input  : - A spectrum to match. This can be given in one of the following
%            forms:
%            (i) A column vector of wavelengths [Ang] of spectral lines in
%            the spectrum to attempt and match with the template.
%            (ii) Two column matrix containing the spectrum.
%            In this case the program will identify the lines in the
%            spectrum.
%            (iii) A FITS file name (wil be read using read_spec.m).
%          - Template to match. This can be one of the following:
%            (i) Column vector of wavelengths [Ang] of spectral lines in
%            the template to attempt and match with the spectrum.
%            (ii) Two column matrix containing the template spectra.
%            In this case the program will identify the lines in the
%            template.
%            (iii) A string containg a template name.
%            For templates availability see spec_get_arc.m.
%          * Arbitrary pairs of ...,keyword,value,...
%            The following keywords are available:
%            'ContrastRMS'  - RMS of contrast in peak finding algorithm.
%                             See find_contrast_peaks.m for details.
%                             Default is 8.
%            'MaxLineShift' - Maximum tolarence between spectrum and template
%                             lines to match. Default is 10.
%            'FitPoly'      - Fit polynomial to matched lines
%                             (pos vs. wavelength).
%                             {'y'|'n'}. Default is 'y'.
%            'PolyFitPar'   - Cell array of additional parameters to pass to
%                             fit genpoly.m for the fotting process.
%                             Default is {'Algo','chol','NormX','y',...
%                             'MaxIter',1,'Method','StdP','Mean','median',...
%                             'Clip',[2 2]};
%            'PlotMS'       - If spectrum is available then will plot the
%                             spectrum, with the template lines, the spectrum
%                             lines and the matched lines. The matched lines
%                             will be divided into two groups.
%                             Lines with a single possible match, and lines
%                             with multiple possible matched (within
%                             threshold). {'y'|'n'}. Default is 'n'.
%            'PlotFit'      - Plot best fit for wavelength calibration,
%                             {'y'|'n'}, default is 'n'.
%            'Deg'          - Degree of polynomial fit. Default is 2.
%            'Interactive'  - Run in interactive mode (the user will be able
%                             to modify the fit. {'y'|'n'}. Default is 'n'.
%            'MarkerSize'   - Marker size in PlotFit. Default is 12.
% Output : - A structure with the polynomail fit information.
%            See fitgenpoly.m for details.
%          - A structure array with the matched lines.
%            Containing the following fields
%            .PossMatch   - number of possible matches for template line.
%            .TempWave    - The wavelength of the line in the template.
%            .SpecWave    - The wavelength of the best matched spectrum line.
%            .IndSM       - Index of template line.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Mar 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

import Util.fit.*


if (nargin==0),
   error('Illegal number of input arguments');
end

%--- set default parameters ---
DefV.ContrastRMS  = 10;
DefV.MaxLineShift = 3;
DefV.FitPoly      = 'y';
DefV.PolyFitPar   = {'Algo','chol','NormX','y','MaxIter',1,'Method','StdP','Mean','median','Clip',[2 2]};
DefV.PlotMS       = 'n';   % plot spectra with matched lines
DefV.PlotFit      = 'n';   
DefV.Deg          = 2;
DefV.Interactive  = 'n';
DefV.MarkerSize   = 12;

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});


%--- Spec ---
if (ischar(Spec)),
   % read from file
   ST   = read_spec(Spec);
   Spec = [ST.Wave, ST.Flux];
else
   if (min(size(Spec))==1),
      % list of wavelength of spectral line is provided.
      SpecLines = Spec;
      Spec      = [];
   else
      % assume spectrum is given - find lines.
      % Identify lines in Spec
      %save try.mat
      [SelectedMaxima,~] = find_contrast_peaks(Spec(:,1),...
						                       Spec(:,2),...
                                               InPar.ContrastRMS);
      SpecLines = SelectedMaxima(:,1);
   end
end

%--- Template ---
if (ischar(Template)),
   % read from file
   AT = spec_get_arc(Template);
   TemplateLines = AT.Lines;
else
   if (min(size(Template))==1),
      % list of wavelength of spectral line is provided.
      TemplateLines = Template;
      Template      = [];
   else
      % assume spectrum is given - find lines.
      % Identify lines in Template
      [SelectedMaxima,~] = find_contrast_peaks(Template(:,1),...
						                       Template(:,2),...
                                               InPar.ContrastRMS);
      TemplateLines = SelectedMaxima(:,1);
   end
end


%--- match lines ---
Iline = 0;
Nselmax = length(TemplateLines);
for Iselmax=1:1:Nselmax,
   [MinLineShift,MinInd] = min(abs(TemplateLines(Iselmax) - SpecLines));
   if (MinLineShift<InPar.MaxLineShift),
      % possible line identified
      Iline = Iline + 1;
      MatchedList(Iline).PossMatch = length(MinInd);   % number of possible matches
      MatchedList(Iline).TempWave  = TemplateLines(Iselmax);         % observed wavelength
      MatchedList(Iline).SpecWave  = SpecLines(MinInd);   % theoretical wavelength
      MatchedList(Iline).IndSM     = Iselmax;
   end
end
if (Iline==0),
    MatchedList = [];
end

% FlagGood  = [MatchedList.PossMatch]==1;
% Resid     = [MatchedList(FlagGood).SpecWave].' - [MatchedList(FlagGood).TempWave].';
% Lambda    = [MatchedList(FlagGood).TempWave].';
% 
% MinLambdaSL = min(SpecLines);
% MaxLambdaSL = max(SpecLines);
% RangeLambdaSL = range(SpecLines);
% 
% Ifit      = find(abs(Resid) < InPar.MaxLineShift);

% if ((RangeLambdaSL.*0.8)>range(Lambda)),
%      % attempt to find more lines by extrapolating relation
%      Par = polyfit(Lambda,Resid,2);
%      Par 
%      

% Par=fitgenpoly(Lambda,Resid,1,2,'Clip',[-2 2],'MaxIter',3);





% Plot matched on top of spectrum
switch lower(InPar.PlotMS)
   case 'y'
       % plot only if spectrum is available
       if (isempty(Spec)),
          if (isempty(Template)),
               % can not plot
               SpPlot = [];
          else
               SpPlot = Template;
          end
       else
	       SpPlot = Spec;
       end
       if (~isempty(SpPlot)),

          figure;
          plot(SpPlot(:,1),SpPlot(:,2),'k-')
          hold on;
          plot(SpecLines,interp1(SpPlot(:,1),SpPlot(:,2),SpecLines),'ro','MarkerSize',InPar.MarkerSize);
          plot(TemplateLines,interp1(SpPlot(:,1),SpPlot(:,2),TemplateLines),'rx','MarkerSize',InPar.MarkerSize);
          Isingle_match = find([MatchedList.PossMatch]==1);
          Imult_match = find([MatchedList.PossMatch]>1);

          plot([MatchedList(Isingle_match).SpecWave],...
               interp1(SpPlot(:,1),SpPlot(:,2),[MatchedList(Isingle_match).SpecWave]),'bx','MarkerSize',InPar.MarkerSize);
          plot([MatchedList(Imult_match).SpecWave],...
               interp1(SpPlot(:,1),SpPlot(:,2),[MatchedList(Imult_match).SpecWave]),'b^','MarkerSize',InPar.MarkerSize);
          legend('Arc','Spectrum lines','Template lines','Matched (single)','Matched (mult)');

          xlabel('Wavelength [A]');
          ylabel('Intensity');
       end

    otherwise
       % do nothing
end


if (numel(MatchedList)<InPar.Deg),
    % no solution can be found
    Poly = [];
else
    % fit polynomial
    
    switch lower(InPar.FitPoly)
     case 'y'
        switch lower(InPar.Interactive)
           case 'n'
               % normalize X-axis
               CenterX = 0; %median([MatchedList.LineWave].');                                           
               X = [MatchedList.TempWave].' - CenterX;
               Y = [MatchedList.SpecWave].'-[MatchedList.TempWave].';
               % fit polynomial
            
               if (length(X)<=InPar.Deg),
                   Poly = [];
               else
                  Poly = fitgenpoly(X,Y,1,InPar.Deg,InPar.PolyFitPar{:});
               end
               
           case 'y'
               % normalize X-axis
               CenterX = 0; %median([MatchedList.LineWave].');       
               X = [MatchedList.TempWave].' - CenterX;
               Y = [MatchedList.SpecWave].'-[MatchedList.TempWave].';
               % plot and fit polynomial
               figure;
               plot(X,Y,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
               H = xlabel('Wavelength [A]');
               set(H,'FontSize',16);
               H = ylabel('Residual (Spec-Temp) [A]');
               set(H,'FontSize',16);

               ResInt = plot_int([],[],...
                                 'fitgenpoly',...
                                 {1,InPar.Deg, InPar.PolyFitPar{:},'Plot','fitonly'},...
                                 'FunParInd',2,...
                                 'FunBehav','i');
               waitfor(gcf,'KeyPressFcn','');

               Poly = ResInt.FunOut{1};
           otherwise
               % do nothing
        end
     case 'n'
        % do not fit lines
        Poly = [];
     otherwise
        error('Unknwon FitPoly option');
    end

    %Fit.LFwave = Fit.XCwave + polyval(Poly.Par,Fit.XCwave);
    %Fit.Poly   = Poly;

    switch lower(InPar.PlotFit)
       case 'y'
          figure;
          plot(X,Y,'ko');
          hold on;
          plot(X(Poly.FlagUse==0),Y(Poly.FlagUse==0),'rx','MarkerSize',12);
          plot(X,polyval(Poly.Par,X),'r-')


       otherwise
          % do nothing
    end
end

