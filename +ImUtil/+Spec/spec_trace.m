function [Trace]=spec_trace(Sim,StartPos,varargin)
%--------------------------------------------------------------------------
% spec_trace function                                               ImSpec
% Description: Trace and extract spectrum from a 2-dimensional image.
% Input  : - The 2-D spectrum image.
%            This can be either a matrix or a file name.
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
%            'PlotTraceIm'- Show the final trace plotted over the 2D
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
%             Extracte             {'no'|'median'|'mean'}. Default is 'no'.
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
%            'UseMask'   - Use image bit mask information
%                          {true|false}. Default is true.
%            'ProblematicBits' - The indices of the bits to flag as bad
%                          bits. Default is [2 9].
%            'FracBadBits' - Maximum allowed frcation of the bad mask
%                          bits within the sliding trace window.
%                          Default is 0.1.
%                          Dispersion positions with fraction of bad bits
%                          larger than this threshold will not be
%                          used in the polynomial fit.
%            'SemiWBadPix' - Semi width in which to count the fraction of
%                          bad pixels. Default is 5.
%            'MasterTrace' - A master trace [Dispersion axis, Spatial axis]
%                          of a trace star. The program is looking for
%                          large deviations between the shape of the
%                          trace and the master trace.
%                          If 'line', then will use a stright line
%                          master trace. Default is 'line'.
%            'MaxDevMT'  - Maximum distance between the master trace
%                          and the trace. Points with distance above this
%                          threshold will not be used in the fit.
%                     xt     Default is 100. When real master trace is
%                          provided it is recomended to decrease this to
%                          5-10.
%            'MinDN'     - Minimum number of counts in the peak of the
%                          trace to use in the fit. Default is 10.
%            'Int'       - Interactive fitting mode {true|false}.
%                          Default is false.
%            'FitType'   - Trace smooth function fitting. Options are:
%                          'poly' - polynomial fit.
%                          'mpoly' - Fit multiple polynomial degrees and
%                                   choose the best in term of maximum
%                                   deviation (defaut).
%                          'MT'   - Master trace fit.
%            'PolyDeg'   - Polynomial degree for the fit
%                          If FitType=poly then default is 3.
%                          If FitType=mpoly then default is (1:1:10).'.
%            'fitgenpolyPar' - Cell array of additional parameters to
%                          pass to fitgenpoly.m.
%                          Default is {'MaxIter',2,'Clip',[2
%                          2],'NormX','no'}.
%            'Bit_TraceDev' - The index of the bit in the bit mask which
%                          to set to true in case that there is a big
%                          deviation (larger than 'MaxOffset') between
%                          the trace and the best fitted
%                          smooth trace. Alternatively, this can be
%                          a function handlde.
%                          Default is @def_bitmask_specpipeline;
%            'MaxOffset' - Maximum offset between the fitted position
%                          and measure positions. Points with larger offset
%                          will be flagged. Default is 1.
%            'MaskType'  - Bit mask type to us, if created.
%                          Default is 'uint16'.
% Output : - Trace structure containing the following fields:
%            .X       - Position along the dispersion direction.
%            .Y       - Traced position along the spatial position axis.
%            .SmY     - Smooth/fitted traced Y position (e.g., polynomial
%                       fit).
%            .Mask    - One D bit mask along the trace, in which the 
%                       'Bit_TraceDev' is set.
%            .StdY    - second moment of the Y position.
%            .H       - Pixel counts at the traced position.
%            .FracBadPix - The fraction of bad pixels within the sliding
%                       window around the traced position. The sliding
%                       window half size is BinDisp*ExtSemiW.
%                       Bad pixels are defined as those for which the
%                       'ProblematicBits' are true.
%            .dYdX    - dY/dX in the trace position. This is useful in
%                       order to identify bad trace.
%            .dYdX5   - dY per 5 pixels in X.
%            .MasterOffset - Offset of the trace relative to the best
%                       fit master trace.
%            .Poly    - A structure containing the best fit polynomial
%                       information. See fitgenpoly.m for details.
%                       Trace.Poly.RMS - contains the best fit rms.
%            .Resid   - Vector of residuals between the observed trace
%                       and fitted trace.
%            .RMS     - The best fit RMS, but including all the pixels
%                       including those that didn't used in the fit.
%            .Chi2    - \chi^2 of best fit (all points).
%            .Flag    - Structure of various flags information.
% See also: spec_trace_int.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Trace,ExtractedSpec]=spec_trace('lred0064.fits',[]);
%          Trace    = spec_trace(Sim(Is),StartPos,'Int',false,'Display','ds9','PlotTraceIm',true,...
%                                                 'MasterTrace',MasterTrace,'MaxDevMT',InPar.MaxDevMT);         
% Reliable: 2
%--------------------------------------------------------------------------

import Util.fit.*
import Util.array.*

Def.StartPos = [];
if (nargin==1)
   StartPos   = Def.StartPos;
end


ImageField  = 'Im';
%HeaderField = 'Header';
%FileField   = 'ImageFileName';
MaskField   = 'Mask';
%BackImField = 'BackIm';
%ErrImField  = 'ErrIm';



% set default values:
% display and interaction
DefV.Display            = 'ds9';   % {'ds9' | 'imtool'}  % need to fix bug with imtool
DefV.DS9Type            = 'new';   % present in new frame
DefV.PlotTraceIm        = false;
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
% masking
DefV.UseMask            = true;
DefV.Mask               = [];   % override the Mask field if available
DefV.ProblematicBits    = [2 9];
DefV.SemiWBadPix        = 5;    % semi width for bad pixel fraction
DefV.MaxFracBadBits     = 0.1;
DefV.MasterTrace        = 'line';  % if empty - do not use
DefV.MaxDevMT           = 100; %5;   % decrease
DefV.MinDN              = 0;  % remove Height less than this value
% bit mask
DefV.Bit_TraceDev       = @def_bitmask_specpipeline;
DefV.MaxOffset          = 1;   % offset between fitted position and measured position threshold
DefV.MaskType           = 'uint16';
% fitting
DefV.Int                = false;
DefV.FitType            = 'mpoly';   % {'poly'|'MT'}
DefV.PolyDeg            = 3;
DefV.fitgenpolyPar      = {'MaxIter',2,'Clip',[2 2],'NormX','no'};


InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Def.MpolyDeg = (1:1:10).';
switch lower(InPar.FitType)
    case 'mpoly'
        if (numel(InPar.PolyDeg)==1)
            InPar.PolyDeg = Def.MpolyDeg;
        end
    otherwise
        % do nothing
end

% get bit mask
InPar.Bit_TraceDev      = get_bitmask_def(InPar.Bit_TraceDev,'Bit_TraceDev');


%-----------------------------------------
%--- Read image and trace the spectrum ---
%-----------------------------------------
[Trace,StartPos,Image]=spec_trace_int(Sim,StartPos,varargin{:});


% set mask
if (InPar.UseMask)
    if (isempty(InPar.Mask))
       if (isfield(Sim,MaskField))
           InPar.Mask = Sim.(MaskField);
       else
           % mask is not available
           InPar.UseMask = false;
           warning('MATLAB:spec_trace:Mask','Mask is not available - do not use mask');
       end  
    end
end
        
% check mask
if (InPar.UseMask)
    FlagBad = maskflag_check(InPar.Mask,sum(2.^(InPar.ProblematicBits-1)),'and');
else
    FlagBad = false(size(Image));
end


% Image parameters
SizeIm    = size(Image);
SizeDisp  = SizeIm(2);
%SizeSpat  = SizeIm(1);
%VecSpat   = (1:1:SizeSpat).';
VecDisp   = (1:1:SizeDisp).';

if (isempty(InPar.GoodRange)
    InPar.GoodRange = [1 SizeDisp];
end

% master trace
if (~isempty(InPar.MasterTrace))
   if (ischar(InPar.MasterTrace))
       switch lower(InPar.MasterTrace)
           case 'line'
               % straight line in dispersion direction
               MasterTrace = [VecDisp, ones(SizeDisp,1)];
           otherwise
               error('Unknown MasterTrace option');
       end
   else
       % interpolate to current grid
       MasterTrace = [VecDisp, interp1(InPar.MasterTrace(:,1),InPar.MasterTrace(:,2),VecDisp,InPar.InterpMethod,'extrap')];
   end
end


% estimate the fraction of bade pixels along the trace
Trace.FracBadPix = zeros(SizeDisp,1);
for Idisp=1:1:SizeDisp
    PosX = Trace.X(Idisp);
    PosY = Trace.Y(Idisp);
    X1 = max(PosX-InPar.BinDisp,1);
    X2 = min(PosX+InPar.BinDisp,SizeDisp); 
    Y1 = round(PosY-InPar.SemiWBadPix);
    Y2 = round(PosY+InPar.SemiWBadPix);
    
    Trace.FracBadPix(Idisp) = sumnd(FlagBad(Y1:Y2,X1:X2))./((X2-X1).*(Y2-Y1));
end
   
Trace.dYdX  = [Inf; diff(Trace.Y)];
Trace.dYdX(isnan(Trace.dYdX)) = Inf;
Trace.dYdX5 = interp1(Trace.X(1:5:end),[Inf; diff(Trace.Y(1:5:end))],Trace.X,InPar.InterpMethod,'extrap');

%plot(Trace.X,Trace.dYdX5,'.')

%[Trace.X, Trace.Y, Trace.H, Trace.dYdX, Trace.dYdX5]
MasterOffset = (Trace.Y - nanmedian(Trace.Y)) - (MasterTrace(:,2) - nanmedian(MasterTrace(:,2)));
Trace.MasterOffset = MasterOffset;

Flag.GoodRange  = Trace.X>=InPar.GoodRange(1) & Trace.X<=InPar.GoodRange(2);
Flag.MT         = abs(MasterOffset)<InPar.MaxDevMT;

Flag.dYdX       = abs(Trace.dYdX)<0.03;
Flag.dYdX5      = abs(Trace.dYdX5)<0.1;
Flag.DN         = Trace.H>InPar.MinDN;
Flag.GoodBitPix = Trace.FracBadPix<InPar.MaxFracBadBits;

FlagAll         = Flag.GoodRange & Flag.MT & Flag.dYdX & Flag.dYdX5 & Flag.DN & Flag.GoodBitPix;




%plot(Trace.X(FlagAll),Trace.Y(FlagAll),'g.')

%save try1.mat

%-------------------------------------
%--- Fit a polynomial to the trace ---
%-------------------------------------
switch lower(InPar.FitType)
    case 'mt'
        Trace.SmY = MasterTrace(:,2) + nanmedian(Trace.Y) - nanmedian(MasterTrace(:,2));
        %plot(Trace.X(FlagAll),Trace.Y(FlagAll),'.')
        %hold on
        %plot(Trace.X,Trace.SmY,'r-')

    case 'poly'
       
        Res = fitgenpoly(Trace.X(FlagAll),Trace.Y(FlagAll),[],...
                         InPar.PolyDeg,InPar.fitgenpolyPar{:});
        %Trace.SmY = polyval(Res.Par,Trace.X);
        Poly = Res;
        Trace.Poly  = Poly;
        [Trace.SmY,~] = polyconf_cov(Poly.Par,Trace.X,Poly.Cov);

    case 'mpoly'
        % fit several polynomials and choose the best
        Npoly = numel(InPar.PolyDeg);
        MaxResid = zeros(Npoly,1).*NaN;
        for Ipoly=1:1:Npoly
            Par   = polyfit(Trace.X(FlagAll),Trace.Y(FlagAll),InPar.PolyDeg(Ipoly));
            PolyY = polyval(Par,Trace.X);
            MaxResid(Ipoly) = max(abs(Trace.Y(FlagAll) - PolyY(FlagAll)));
        end
        %plot(MaxResid)
        [MinMaxResid,MinInd] = min(MaxResid);
        InPar.PolyDeg = InPar.PolyDeg(Ipoly);
        Par   = polyfit(Trace.X(FlagAll),Trace.Y(FlagAll),InPar.PolyDeg);
        
        
        %Res = fitgenpoly(Trace.X(FlagAll),Trace.Y(FlagAll),[],...
        %                 InPar.PolyDeg,InPar.fitgenpolyPar{:});
        %Trace.SmY = polyval(Res.Par,Trace.X);
        Trace.Poly.Par  = Par;
        Trace.SmY       = polyval(Par,Trace.X);
        Trace.Resid     = Trace.Y - Trace.SmY;

        Trace.Poly.RMS  = nanstd(Trace.Resid(FlagAll));  

        %[Trace.SmY,~] = polyconf_cov(Poly.Par,Trace.X,Poly.Cov);

        
    otherwise
        error('Unknown FitType option');
end
Trace.Resid    = Trace.Y - Trace.SmY;

Flag.BigOffset = abs(Trace.Resid)>InPar.MaxOffset;
Trace.RMS      = nanstd(Trace.Resid);

if (InPar.Int)
    % Interactive rejection of points
    %MaxIter    = 1;
    %PlotOption = 1;
    %[XY2,XY1,XY0,Par,Res,IU,INU,GraphicH] = plot_rm(Trace.Disp,Trace.Spat,'b.','polyfit_sc',5,InPar.Deg,InPar.SigClip,MaxIter,PlotOption);

    plot(Trace.X(FlagAll),Trace.Y(FlagAll),'o');
    H = xlabel('Dispersion axis [pix]');
    set(H,'FontSize',16);
    H = ylabel('Spatial axis [pix]');
    set(H,'FontSize',16);
    [Res] = plot_int([],[],...
                     @Util.fit.fitgenpoly,...
                     {[],InPar.PolyDeg,InPar.fitgenpolyPar{:},'Plot','fitonly'},...
                     'FunParInd',InPar.PolyDeg,...
                     'FunBehav','i');

    %waitfor(gcf,'KeyPressFcn','');
    Poly = Res.FunOut{1};
    Trace.Poly  = Poly;
    [Trace.SmY,~]=polyconf_cov(Poly.Par,Trace.X,Poly.Cov);

end

Trace.Flag = Flag;

% plot trace on image
if (InPar.PlotTraceIm)
    switch lower(InPar.Display)
        case 'ds9'
            ds9_disp(Sim.(ImageField));
            ds9_plotregion(Trace.X,Trace.SmY,'Type','line','Size',[Trace.X(1:end-1), Trace.SmY(1:end-1), Trace.X(2:end), Trace.SmY(2:end)]);
            
        case 'imtool'
            imtool(Sim.(ImageField));
            plot(Trace.X,Trace.SmY,'r-');
            
        otherwise
            error('Unknown Display option');
    end
end
    
% wtite bit mask
if (~isempty(InPar.Bit_TraceDev))
   Trace.Mask = maskflag_set([],InPar.MaskType, InPar.Bit_TraceDev, Trace.Resid>InPar.MaxOffset);
end



    
