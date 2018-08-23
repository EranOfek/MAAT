function [TracePeaks,Spec]=spec_select_tracepoints(FileName,varargin)
%------------------------------------------------------------------------------
% spec_select_tracepoints function                                      ImSpec
% Description: Given a spectrum in a matrix or FITS file format, collapse
%              the spectrum in the dispersion axis and search for all
%              the traces in the spectrum. The output of this function
%              is the input for spec_trace.m
% Input  : - A matrix containing the spectrum, or a FITS file name (string).
%          * Arbitrary number of pairs of ...,keyword,value,... arguments.
%            The following keywords are available:
%            'DispAxis'    - dispesion axis {'x','y'}. Default is 'x'.
%                            Alternatively, this can be a keyword header,
%                            that contains either 'x' or 'y'.
%                            If empty matrix, then the program will guess
%                            the dispersion axis by assuming it is
%                            the longest axis of the image.
%            'AlgoPeaks'   - Algorithm to identify traces:
%                            'single' - select peaks along a single cut
%                                       in the position direction (default).
%            'CutSemiWidth'- Semi width of cut. This is the half-width of
%                            the cut along the position direction by which
%                            to median the data before peaks are searched.
%                            Default is 100.
%            'AlgoCollapse'- The collapse algorith
%                            {'mean','median'}. Default is 'median'.
%            'CutPosition' - The position along the dispersion axis of the
%                            center of the cut. This can be an integer
%                            specifying a position along the dispersion
%                            axis, or a string 'middle'
%                            (middle of the spectrum) - default.
%            'ContrastRMS' - Threshold contrast to use when searching for 
%                            peaks in the collapsed region along the
%                            position ais. This is given in units of
%                            the rms (default is 10).
%                            This parameter refers to the rms of the peak
%                            in the collapsed data, rather than the rms in
%                            a single pixel.
%            'AlgoPrms     - Algorithm by which to identify peaks
%                            (see find_contrast_peaks.m). Default is 'win'.
%            'StdPar'      - Additional parameters to pass to
%                            find_contrast_peaks.m. Default is
%                            {30,[],0.68}.
%            'FitsPar'     - A cell array containing additional parameters
%                            to pass to fitsread.m (e.g., {'primary',1}).
%                            Default is {}.
%            'Manual'      - Selecting peaks manualy, or semi manualy.
%                            This can be one of the following options:
%                            'auto' - automatic peaks selection (default).
%                            'plot' - will disply the collapsed cut along
%                                     the position direction, and the user
%                                     will be prompted to select or reject
%                                     peaks. 
%                            'ds9'  - will open ds9 (if not opened) and
%                                     display the image with the possible
%                                     peaks marked. The user will be
%                                     prompted to select traces
%                                     The program will simply return the
%                                     user selection without refinment.
%            'Plot'        - {'y'|'n'} Open a new figure and presenting
%                            the collapsed cut along the position direction
%                            with the selected peaks marked. Default is 'n'.
% Output : - A structure containing information about the possible peaks
%            that can be traced in the spectrum. The following fields
%            are available.
%            .PeakPos       - The position of the peak along the position
%                             direction.
%            .PeakDispPos   - The position of the peak along the disperssion
%                             position (i.e., the cut center).
%            .PeakHeight    - The height of the peak.
%            .PeakSN        - The peak signal/noise.
%          - A matrix containing the spectrum.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Dec 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [TracePeaks,Spec]=spec_select_tracepoints('lred0131.fits');
%          [TracePeaks,Spec]=spec_select_tracepoints('lred0131.fits','Manual','ds9');
%------------------------------------------------------------------------------

DefV.DispAxis     = 'x';
DefV.AlgoPeaks    = 'single';
DefV.CutSemiWidth = 100;
DefV.AlgoCollapse = 'median';
DefV.CutPosition  = 'middle';
DefV.ContrastRMS  = 10;
DefV.AlgoPrms     = 'win';
DefV.StdPar       = {30,[],0.68};
DefV.FitsPar      = {};
DefV.Manual       = 'auto';
DefV.Plot         = 'n';

InPar = set_varargin_keyval(DefV,'y','use',varargin{:});


%--- Read Spectrum from file ---
if (ischar(FileName)),
   FileType = 'fits';
   Spec = fitsread(FileName,InPar.FitsPar{:});
else
   FileType = 'matrix';
   Spec = FileName;
   clear FileName;
end


%--- Init dispersion axis direction ---
% find or define the dispersion axis on image
switch lower(InPar.DispAxis)
 case 'x'
    ColXY = [1 2];
 case 'y'
    ColXY = [2 1];
 otherwise
    % look for keyword header with this name
    switch lower(FileType)
     case 'fits'
        [KeywordVal]=get_fits_keyword(FileName,InPar.DispAxis);
        if (isnan(KeywordVal{1})),
	       error(sprintf('Header keyword %s doesnt exist',InPar.DispAxis));
        else
	       InPar.DispAxis = KeywordVal{1};
           if (isnumeric(InPar.DispAxis)),
	          % in this case DispAxis can be 1 (for x), or 2 (for y).
              StrVec = {'x','y'};
              InPar.DispAxis = StrVec{InPar.DispAxis};
           end
        end
     case 'matrix'
        error('DispAxis is given as header keyword when input is a matrix');
    end
end

%--- Image properties ---
% image size
[Ny,Nx] = size(Spec);



%--- peaks selection algorithm ---
% look for peaks along the position direction
% after collapsing a section of the image
switch lower(InPar.AlgoPeaks)
 case 'single'

    % find location in which to cut the image
    if (ischar(InPar.CutPosition)),
       switch lower(InPar.CutPosition)
        case 'middle'
           switch lower(InPar.DispAxis)
            case 'x'
               InPar.CutPosition = floor(Nx./2);
            case 'y'
               InPar.CutPosition = floor(Ny./2);
            otherwise
	       error('Unknwon DispAxis option');
           end
        otherwise
	   error('Unknwon CutPosition option');
       end
    end

    % cut the image and verify that the cut is not out of image bounds
    switch lower(InPar.DispAxis)
     case 'x'
        Start   = max([1;InPar.CutPosition-InPar.CutSemiWidth]);
        End     = min([Nx;InPar.CutPosition+InPar.CutSemiWidth]);
        CutSpec = Spec(:,Start:End);
        DimCollapse = 2;
     case 'y'
        Start   = max([1;InPar.CutPosition-InPar.CutSemiWidth]);
        End     = min([Ny;InPar.CutPosition+InPar.CutSemiWidth]);
        CutSpec = Spec(Start:End,:);
        DimCollapse = 1;
     otherwise
        error('Unknwon DispAxis option');
    end

    % collapse the image
    switch lower(InPar.AlgoCollapse)
     case 'median'
        Collapse = nanmedian(CutSpec,DimCollapse);
     case 'mean'
        Collapse = nanmean(CutSpec,DimCollapse);
     otherwise
        error('Unknwon AlgoCollapse option');
    end

    % search for peaks in collapse vector
    PosAxis = (1:1:length(Collapse))';
    [SelectedMaxima,Cont_RMS] = find_contrast_peaks(PosAxis,...
						    Collapse,...
						    InPar.ContrastRMS,...
                            InPar.AlgoPrms,...
                            InPar.StdPar);

    % sort out the peaks and their properties
    TracePeaks.PeakPos     = SelectedMaxima(:,1);
    TracePeaks.PeakHeight  = SelectedMaxima(:,2);
    TracePeaks.PeakDispPos = InPar.CutPosition.*ones(size(SelectedMaxima(:,1)));
    TracePeaks.PeakSN = Cont_RMS;


 otherwise
    error('Unknown AlgoPeaks option');
end


%--- Manual operation/peaks selection ---
% allow the user to select or remove peaks
% using either a plot of the collapse image or displaying
% the full image in ds9.
switch lower(InPar.Manual)
 case 'auto'
    % do nothing - automatic procedure completed
 case 'ds9'
    % start ds9 (if not yet started)
    % this command will work only in Linux/Unix based systems
    ds9_start;

    % if FITS image is not available then save a FITS image
    switch lower(FileType)
     case 'matrix'
        FileName = 'tmp.spec_select_tracepoints.fits';
        fitswrite(Spec,FileName);
     otherwise
        % do nothing
    end

    % display the image
    ds9_disp(FileName);

    % mark possible peaks in ds9 display.
    Cat = [TracePeaks.PeakDispPos, TracePeaks.PeakPos];
    [Cat,Col]=ds9_dispsex(Cat,ColXY);

    clear TracePeaks;
    % prompt the user to select peaks
    Select = 1;
    TraceInd = 0;
    while (Select==1)
       TraceInd = TraceInd + 1;
       fprintf('Select trace number %d\n',TraceInd);
       %fprintf('r - remove; n - stop; y - select\n');
       Suser = input('Use the mouse (left click) to select a trace [y]/n : ','s');
       switch lower(Suser)
        case 'n'
           Select = 0;
        otherwise
           [CooX,CooY,Value]=ds9_getcoo(1,'image');
           DS9XY = [CooX, CooY];

           TracePeaks.PeakPos(TraceInd)     = DS9XY(ColXY(2));
           TracePeaks.PeakHeight(TraceInd)  = Value;
           TracePeaks.PeakDispPos(TraceInd) = DS9XY(ColXY(1));
           TracePeaks.PeakSN(TraceInd)      = NaN;
       end
    end
 case 'plot'
    % plot the Collapse vector and select peaks
    stairs((1:1:length(Collapse)).',Collapse,'k-');
    % plot possible peaks
    hold on;
    %plot(TracePeaks.PeakPos,TracePeaks.PeakHeight,'rv','MarkerSize',12);

    fprintf('Select or remove peaks\n');
    [Res]=plot_interactive(TracePeaks.PeakPos,...
                           TracePeaks.PeakHeight,...
		                   {'rv','MarkerSize',12},{},...
		                   'Position [pix]','Height [DN]');

    TracePeaks.PeakPos     = TracePeaks.PeakPos(Res.Ind);
    TracePeaks.PeakHeight  = TracePeaks.PeakHeight(Res.Ind);
    TracePeaks.PeakDispPos = TracePeaks.PeakDispPos(Res.Ind);
    TracePeaks.PeakSN      = TracePeaks.PeakSN(Res.Ind);
 otherwise
    error('unknown Manual option');
end 


%--- Plot --- 
% plot collapse and peaks
switch lower(InPar.Plot)
 case 'y'
    stairs((1:1:length(Collapse)).',Collapse,'k-');
    % plot possible peaks
    hold on;
    plot(TracePeaks.PeakPos,TracePeaks.PeakHeight,'rv','MarkerSize',12);
 case 'n'
    % do nothing
 otherwise
    error('Unknown Plot option');
end


    




