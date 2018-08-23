function [Trace,StartPos,Image]=spec_trace_int(Sim,StartPos,varargin)
%--------------------------------------------------------------------------
% spec_trace_int function                                           ImSpec
% Description: Trace a spectrum in a single image by tracing the peak of 
%              the intensity pixel by pixel in the dispersion direction.
% Input  : - A single image in which to trace the spectrum.
%            The following inputs are possible:
%            (1) Cell array of image name in string format.
%            (2) String containing image name (e.g., 'lref0127.fits').
%            (3) Structure of image (SIM).
%                The image should be stored in the 'Im' field.
%                This may contains also mask image (in the 'Mask' field),
%                and an error image (in the 'ErrIm' field).
%            (4) Cell array containing a matrix.
%            (5) A file contains a single image (e.g., '@list').
%          - [X, Y] points in which to start the spectrum tracing.
%            if empty matrix (i.e., []), then use the 'Display' to
%             choose trace starting position interactively.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Display'   - Display image using {'ds9'|'imtool'}, default
%                          is 'ds9'.
%            'DS9Type'   - In which ds9 frame the image is presented
%                          {'new'|'first'|'last'|'next'|'prev'}.
%                          Default is 'new'.
%                          spectrum {true|false}. Default is false.
%            'ExtSemiW'  - Semi width along the spatial dimension in
%                          which to look and refine for the trace position.
%                          Default is 50.
%            'DispDir '  - Dispesrion axis {'x'|'y'}. Default is 'x'.
%            'BinDisp'   - Number of pixels (in half-bin) to bin in the
%                          dispersion direction (for improving S/N),
%                          default is 20.
%            'StepDisp'  - Reclocating trace position every 'StepDisp'
%                          pixels. Default is 5.
%            'GoodRange' - range along to dispersion axis [Min Max],
%                          to use as a good region.
%                          Default is empty. If empty use all.
%            'TraceCollapse' - Method to use in order to collapse the
%                          spectrum (in the bin 'BinDisp' along the
%                          dispersion direction). Options are:
%                          {'median'|'mean'}. Default is 'median'.
%            'LocalBack' - Local background removal method
%                          {'no'|'median'|'mean'}. Default is 'no'.
%            'PeakMethod'- Trace center (in spatial position) finding
%                          method. Options are:
%                          'wmean' - Mean weighted by intesnity.
%                          'wmeanr'- Mean weighted by intesnity with
%                                    regularization (forcing the trace
%                                    center to be not too far away from
%                                    the position found at the previous
%                                    step). Default.
%            'FirstTimeReg'- Reglarize on the first iteration
%                          {true|false}. Default is false.
%            'InterpMethod'- Interpolation method. Default is 'linear'.
%            'SigmaReg'  - Sigma of the regularization Gaussian.
%                          Default is 5.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
%            image2sim.m
% Output : - A structure containing the trace information. The following
%            fields are available:
%            .X - X (along the dispersion position).
%            .Y - Y position of the trace.
%            .StdY - An indicator of the Y position error.
%            .H - The height at the trace peak.
%          - [X,Y] start position for trace.
%          - Background subtracted image in which the spectrum was trace.
% See also: spec_trace.m
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Trace,StartPos,Image]=spec_trace_int('red0036.fits',[]);
% Reliable: 2
%--------------------------------------------------------------------------


ImageField  = 'Im';
%HeaderField = 'Header';
%FileField   = 'ImageFileName';
%MaskField   = 'Mask';
%BackImField = 'BackIm';
%ErrImField  = 'ErrIm';



% set default values:
% display and interaction
DefV.Display            = 'ds9';   % {'ds9' | 'imtool'}  % need to fix bug with imtool
DefV.DS9Type            = 'new';   % present in new frame
%DefV.PlotTraceIm        = false;
% traceing
DefV.ExtSemiW           = 50;
DefV.DispDir            = 'x';
DefV.BinDisp            = 20;
DefV.StepDisp           = 5;
DefV.GoodRange          = [];
DefV.TraceCollapse      = 'median';
DefV.PeakMethod         = 'wmeanr';
DefV.FirstTimeReg       = false;
DefV.InterpMethod       = 'linear';
DefV.LocalBack          = 'no';
DefV.SigmaReg           = 4;

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});


% read the image
Sim = image2sim(Sim,varargin{:});



% select starting position manually
if (isempty(StartPos)),
    % select trace intial position manually
    switch lower(InPar.Display)
        case 'ds9'
            ds9_disp(Sim.(ImageField),InPar.DS9Type);
            disp('Select a trace starting point in ds9 display');
            [X,Y] = ds9_getcoo(1,'image');
        case 'imtool'
            imtool(Sim.(ImageField));
            disp('Select a trace starting point in imtool display');
            [X,Y] = ginput(1);
        otherwise
            error('Unknown Display option');
    end
    StartPos = [X, Y];     
end
    
% store image
Image = Sim(1).(ImageField);

% rotate image
switch lower(InPar.DispDir)
 case 'x'
    % do nothing
 case 'y'
    % rotate
    Image    = Image.';
    StartPos = [StartPos(2), StartPos(1)];  % [X, Y]
 otherwise
    error('Unknown DispAxis option');
end

StartPosX = round(StartPos(1));
StartPosY = round(StartPos(2));



% subtract background
[~,Image] = image_background(Image,varargin{:});



% Image parameters
SizeIm    = size(Image);
SizeDisp  = SizeIm(2);
%SizeSpat  = SizeIm(1);
%VecSpat   = (1:1:SizeSpat).';
VecDisp   = (1:1:SizeDisp).';

if (isempty(InPar.GoodRange)),
    InPar.GoodRange = [1 SizeDisp];
end


 
%--------------------------
%--- Trace the spectrum ---
%--------------------------

% initalize trace structure
Trace.X      = VecDisp;                       % X position
Trace.Y      = zeros(size(VecDisp)).*NaN;     % best Y position
Trace.StdY   = zeros(size(VecDisp)).*NaN;     % error in Y position
Trace.H      = zeros(size(VecDisp)).*NaN;     % Trace height
%Trace.FracBadPix = zeros(size(VecDisp)).*NaN;

%PeakSpatPos     = StartPosSpat;
LastPeakPosY    = round(StartPosY);

% Trace the rest of the spectrum spectrum
% go backward

Ind = 0;
for PosX=StartPosX:-InPar.StepDisp:InPar.GoodRange(1),
    Ind = Ind + 1;    
        
    X1 = max(round(PosX-InPar.BinDisp),1);
    X2 = min(round(PosX+InPar.BinDisp),SizeDisp);
    Y1 = round(LastPeakPosY-InPar.ExtSemiW);
    Y2 = round(LastPeakPosY+InPar.ExtSemiW);

    %[PosX, X1,X2, Y1, Y2]
    
    [MatX,MatY] = meshgrid(X1:X2,Y1:Y2);
  
    SubImage  = Util.array.nangetind(Image,MatY,MatX);

    % TraceCollapse - e.g., @nanmedian
    switch lower(InPar.TraceCollapse)
        case 'median'
            Cut = nanmedian(SubImage,2);
        case 'mean'
            Cut = nanmean(SubImage,2);
        otherwise
            error('Unknown Trace collapse option');
    end
   
    switch lower(InPar.LocalBack)
        case 'median'
            Cut = Cut - nanmedian(Cut);
        case 'mean'
            Cut = Cut - nanmean(Cut);
        case {'no','none'}
            % do nothing
        otherwise
            error('Unknown Trace collapse option');
    end
    
    VecY = (Y1:1:Y2).';
    switch lower(InPar.PeakMethod)
        case 'wmean'
            % weighted mean
            PeakY = nansum(Cut.*VecY)./nansum(Cut);
            StdY  = sqrt(nansum(Cut.*VecY.^2)./nansum(Cut) - PeakY.^2);
        case 'wmeanr'
            % regularized weighted mean
            if (Ind==1 && ~InPar.FirstTimeReg),
                PeakY = nansum(Cut.*VecY)./nansum(Cut);
                StdY  = sqrt(nansum(Cut.*VecY.^2)./nansum(Cut) - PeakY.^2);
            else
                Regularization = exp(-(VecY-LastPeakPosY).^2./(2.*InPar.SigmaReg.^2));
                PeakY = nansum(Cut.*Regularization.*VecY)./nansum(Cut.*Regularization);
                StdY  = sqrt(sum(Cut.*Regularization.*VecY.^2)./nansum(Cut.*Regularization) - PeakY.^2);
            end
        otherwise
            error('Unknown PeakMethod option');
    end
    if (isnan(PeakY)),
        % set to last peak value
        PeakY = LastPeakPosY;
    else
        LastPeakPosY = PeakY;
    end
    %[PosX, PeakY]

    Trace.Y(PosX)    = PeakY;
    Trace.StdY(PosX) = StdY;
    Trace.H(PosX)    = interp1(VecY,Cut,PeakY,InPar.InterpMethod);
    
end

LastPeakPosY    = round(StartPosY);

Ind = 0;
for PosX=StartPosX:+InPar.StepDisp:InPar.GoodRange(2),
    Ind = Ind + 1;        

    %PosX
    X1 = max(round(PosX-InPar.BinDisp),1);
    X2 = min(round(PosX+InPar.BinDisp),SizeDisp); 
    Y1 = round(LastPeakPosY-InPar.ExtSemiW);
    Y2 = round(LastPeakPosY+InPar.ExtSemiW);
    
    %[X1, X2, Y1, Y2]
    [MatX,MatY] = meshgrid(X1:X2,Y1:Y2);
   
    SubImage  = Util.array.nangetind(Image,MatY,MatX);

    % TraceCollapse - e.g., @nanmedian
    switch lower(InPar.TraceCollapse)
        case 'median'
            Cut = nanmedian(SubImage,2);
        case 'mean'
            Cut = nanmean(SubImage,2);
        otherwise
            error('Unknown Trace collapse option');
    end
    
    switch lower(InPar.LocalBack)
        case 'median'
            Cut = Cut - nanmedian(Cut);
        case 'mean'
            Cut = Cut - nanmean(Cut);
          case {'no','none'}
            % do nothing    
        otherwise
            error('Unknown Trace collapse option');
    end
    
    VecY = (Y1:1:Y2).';
    switch lower(InPar.PeakMethod)
        case 'wmean'
            PeakY = nansum(Cut.*VecY)./nansum(Cut);
            StdY  = sqrt(nansum(Cut.*VecY.^2)./nansum(Cut) - PeakY.^2);
        case 'wmeanr'
            % regularized weighted mean
            if (Ind==1 && ~InPar.FirstTimeReg),
                PeakY = nansum(Cut.*VecY)./nansum(Cut);
                StdY  = sqrt(nansum(Cut.*VecY.^2)./nansum(Cut) - PeakY.^2);
            else
                Regularization = exp(-(VecY-LastPeakPosY).^2./(2.*InPar.SigmaReg.^2));
                PeakY = nansum(Cut.*Regularization.*VecY)./nansum(Cut.*Regularization);
                StdY  = sqrt(nansum(Cut.*Regularization.*VecY.^2)./nansum(Cut.*Regularization) - PeakY.^2);
            end
            
        otherwise
            error('Unknown PeakMethod option');
    end
    if (isnan(PeakY)),
        % set to last peak value
        PeakY = LastPeakPosY;
    else
        LastPeakPosY = PeakY;
    end
    
    Trace.Y(PosX)    = PeakY;
    Trace.StdY(PosX) = StdY;
    Trace.H(PosX)    = interp1(VecY,Cut,PeakY,InPar.InterpMethod);
    
end


% interpolate missing values
Trace.Y    = interp1_nan(Trace.X,Trace.Y,InPar.InterpMethod,'extrap');
Trace.H    = interp1_nan(Trace.X,Trace.H,InPar.InterpMethod,'extrap');
Trace.StdY = interp1_nan(Trace.X,Trace.StdY,InPar.InterpMethod,'extrap');

