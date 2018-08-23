function [Trace,ExtractedSpec]=spec_trace1(Image,StartPos,varargin)
%--------------------------------------------------------------------------
% spec_trace function                                               ImSpec
% Description: Trace and extract spectrum from a 2-dimensional image.
% Input  : - The 2-D spectrum image.
%            This can be either a matrix or a file name.
%          - [X, Y] points in which to start the spectrum tracing.
%            if empty matrix (i.e., []), then choose interactively.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'DispAxis'  - Dispesrion axis {'x'|'y'}, default is 'x'.
%            'Disp'      - Display image using {'ds9'}, default
%                          is 'ds9' (in this case ds9 should be open).
%            'Int'       - Interactive mode {true|false}. If 'n', then try
%                          to use non-interactive mode whenever possible.
%                          Default is 'y'.
%            'FitsPars'  - Additional parameters to pass to the fitsread.m
%                          program. Default is {}.
%            'BinDisp'   - Number of pixels (in half-bin) to bin in the
%                          dispersion direction (for improving S/N),
%                          default is 50.
%            'StepWave'  - Steps (in pixels) in dispersion direction
%                          in which to trace the spectrum. Default is 5.
%            'MasterTrace' - Master trace is preferably a trace of a
%                          trace object (e.g., a bright std star).
%                          This can be used to look for inconsistecies,
%                          or to force the trace along the master trace.
%            'Kernel'    - If not empty, then run a smoothing kernel on
%                          before looking for the peak.
%                          e.g., fun_gauss([1 5 5],(1:1:10).').
%                          Default is empty.
%            'PeakMethod'- Method to loof for the trace peak.
%                          See find_peak_center.m for options.
%                          Default is 'strn'.
%            'PeakSWidth'- Semi width of region in which to search for
%                          the trace peak.
%            'InterpMethod'- Interpolation method. See interp1.m for
%                          options. Default is 'linear'.
%            'FitType'   - Trace fitting method. Options are:
%                          'mt' - fit the master trace.
%                          'poly' - polynomial fitting.
%            'PolyDeg'   - Polynomial degree to fit. Default is 2.
%            'fitgenpolyPar' - Cell array of parameters to pass to
%                          fitgenpoly.m.
%                          Default is {'MaxIter',2,'Clip',[2
%                          2],'NormX','no'}.
%            'MaskType'  - Mask type. Default is 'uint16'.
%            'Bit_BigOffsetMT' - Bit index containing the "BigOffsetMT"
%                          flag. This flag indicates if there is a
%                          deviation in spatial direction between 
%                          the master trace and the measured trace is
%                          larger than 'MaxDevMT'. Default is 1.
%            'Bit_NearEdge' - Bit index containing the "NearEdge" flag.
%                          Indicating if the pixel is close to edges
%                          within the 'RemoveEdges' distance.
%                          Default is 2.
%            'Bit_BigErr' - Bit index containing the "BigErr" flag.
%                          Indicating if the estimated Poisson error
%                          on the peak of the trace is larger than
%                          the relative error indicaed in 'MaxRelErr'.
%                          Default is 3.
%            'Bit_BigOffsetFit' - Bit index containing the "BigOffsetFit"
%                          flag. Indicates that the offset between fitted
%                          position and measured position is larger than
%                          'MaxOffset'. Default is 4.
%            'MaxDevMT'  - Maximum deviation in spatial direction between 
%                          the master trace and the measured trace.
%                          Pixels with larger deviation will not be used
%                          in the trace fitting process and will be
%                          marked in the mask(bit Bit_BigOffsetMT).
%                          Default is 5 pixels.
%            'RemoveEdges' - Number of pixels from the dispersion edges,
%                          not to use in the various fits. Default is 20.
%            'MaxRelErr' - Maximum relative error to use in fitting.
%                          Default is 0.2.
%            'MaxOffset' - Maximum offset between fitted position and
%                          measured position above to set the
%                          "Bit_BigOffsetFit" bit. Default is 1.
%            'ExtSemiW'  - Semi width of extracted spectrum
%                          (see spec_extract_2d.m for details).
%                          Default is 50.
%            'BinSpatSize' - Binning in spatial direction
%                          (see spec_extract_2d.m for details).
%                          Default is 1.
%            'ExMethod'  - Extraction method {'y' | 'vert'}.
%                          (see spec_extract_2d.m for details).
%                          Default is 'vert'.
% Output : - Trace structure containing the following fields:
%            .Disp    - Position along the dispersion direction.
%            .Spat    - Spatial position of the best traced spectrum,
%                       (corresponding to .Disp).
%            .SpatS   - A smooth version of .Spat - smoothed using
%                       medfilt1 with 'BinDisp' window size.
%            .SpatP   - Spatial position of the best traced spectrum,
%                       after polynomial fitting (corresponding to .Disp).
%            .Par     - Best fit polynomial coefficients.
%            .RMS     - Best fit RMS in the spatial direction.
%          - Matrix of extracted 2D spectrum, each row corresponds to
%            X in BestTrace.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Trace,ExtractedSpec]=spec_trace(Files(I).Image,StartPos,IntrumentSetup(Iinst).Prog.spec_trace.Par{:});
%          [Trace,ExtractedSpec]=spec_trace('lred0064.fits',[]);
% Reliable: 
%--------------------------------------------------------------------------
import Util.array.*
import Util.fit.*

Def.StartPos = [];
if (nargin==1)
   StartPos   = Def.StartPos;
end

ImageField  = 'Im';

% set default values:
DefV.DispAxis           = 'x';
DefV.Disp               = 'ds9';
DefV.Int                = true;
DefV.BinDisp            = 20;
DefV.StepWave           = 5;
DefV.MasterTrace        = [(1:1:4096).',ones(4096,1)];
DefV.Kernel             = []; %fun_gauss([1 5 3],(1:1:10).');
DefV.PeakMethod         = 'strn'; %'strn';
DefV.PeakSWidth         = 5;
%DefV.PeakType           = 'nearest';
DefV.InterpMethod       = 'linear';
DefV.FitType            = 'poly';   % {'poly'|'MT'}
DefV.PolyDeg            = 3;
%DefV.PolyPar            = {};
DefV.fitgenpolyPar      = {'MaxIter',2,'Clip',[2 2],'NormX','no'};
DefV.FitsPars           = {};  % fitsread.m additional parameters
%DefV.PolyPar            = {'Algo','chol','NormX','y','MaxIter',2,'Method','StdP','Mean','median','Clip',[2 2]};
DefV.MaskType           = 'uint16';
DefV.Bit_BigOffsetMT    = 1;
DefV.Bit_NearEdge       = 2;
DefV.Bit_BigErr         = 3;
DefV.Bit_BigOffsetFit   = 4;  % big offset between fitted position and measured position
DefV.MaxDevMT           = 100; %5;   % decrease
DefV.RemoveEdges        = 30;  % number of pixels to remove from spectrum edge during the polynomial fit
DefV.MaxRelErr          = 0.2; % Maximum relative error to use in fitting 
DefV.MaxOffset          = 1;   % offset between fitted position and measured position to mark "Bit_BigOffsetFit".
DefV.MinDN              = 10;  % remove Height less than this value
% spectral extraction parameters
DefV.ExtSemiW           = 50;
DefV.BinSpatSize        = 1;
%DefV.InterpMethod       = 'linear';
DefV.ExMethod           = 'vert';     % {'y' | 'vert'}

InPar = InArg.populate_keyval(DefV,varargin,mfilename);



%--- read image ---
%Sim = image2sim(Image);
%ImMat = Sim.(ImageField);

if (ischar(Image)==1)
   ImMat    = fitsread(Image,InPar.FitsPars{:});
   IsTmp    = 0;
else
   ImMat    = Image;
   % write temporay FITS image
   Image    = 'tmp.spec_trace.fits';
   fitswrite(ImMat,Image);
   IsTmp    = 1;
end


if (isempty(StartPos))
   % select StartPos interactively

   switch lower(InPar.Disp)
    case 'ds9'
       ds9_disp(Image,'new');
       disp('Select a trace starting point in ds9 display');
       [X,Y] = ds9_getcoo(1,'image');
       StartPos = [X, Y];     
    
    otherwise
       error('Unknown Disp option');
   end   
end


% rotate image
switch lower(InPar.DispAxis)
 case 'x'
    % do nothing
 case 'y'
    % rotate
    ImMat    = ImMat.';
    StartPos = [StartPos(2), StartPos(1)];
 otherwise
    error('Unknown DispAxis option');
end

% round start position
StartPosDisp = round(StartPos(1));
StartPosSpat = round(StartPos(2));

% Image parameters
SizeIm    = size(ImMat);
SizeDisp  = SizeIm(2);
SizeSpat  = SizeIm(1);
VecSpat   = (1:1:SizeSpat).';
VecDisp   = (1:1:SizeDisp).';

%--------------------------
%--- Trace the spectrum ---
%--------------------------

% initalize trace structure
Trace.Disp   = VecDisp;                       % X position
Trace.Spat   = zeros(size(VecDisp)).*NaN;     % Y position
Trace.SpatP  = zeros(size(VecDisp)).*NaN;     % fitted Y position
Trace.Height = zeros(size(VecDisp)).*NaN;     % Trace height

% find the spectrum peak
% cut spectrum in the direction perpendicular to the dispersion axis
%[MinInd,MaxInd] = Util.array.check_range(SizeDisp,StartPosDisp-InPar.BinDisp,StartPosDisp+InPar.BinDisp);

% cut along the spatial axis
%CutSpat  = median(ImMat(:,MinInd:MaxInd),2);

% find peak near StartPos  

PeakSpatPos     = StartPosSpat;
LastPeakSpatPos = PeakSpatPos;

% Trace the rest of the spectrum spectrum
Istart = StartPosDisp;
% go backward
for Ipos=Istart:-InPar.StepWave:1
   %Ipos
   [MinInd,MaxInd] = Util.array.check_range(SizeDisp,Ipos-InPar.BinDisp,Ipos+InPar.BinDisp);
   CutSpat  = nanmedian(ImMat(:,MinInd:MaxInd),2);
   
   RegionInd = (floor(LastPeakSpatPos-InPar.PeakSWidth):1:ceil(LastPeakSpatPos+InPar.PeakSWidth)).';
   
   if (~isempty(InPar.Kernel))
      CutSpat  = abs(ifft(fft(CutSpat).*fft(InPar.Kernel,length(CutSpat))));
   end
   
   [PeakSpatPos,PeakH]=find_peak_center(VecSpat,CutSpat,'Method',InPar.PeakMethod,'Ind',RegionInd);
   if (~isnan(PeakSpatPos))
       LastPeakSpatPos = PeakSpatPos;
   end
   Trace.Spat(Ipos) = PeakSpatPos;
   Trace.Height(Ipos) = PeakH;
end

LastPeakSpatPos = StartPosSpat;

% go forward
for Ipos=Istart:InPar.StepWave:SizeDisp
   %Ipos
   [MinInd,MaxInd] = Util.array.check_range(SizeDisp,Ipos-InPar.BinDisp,Ipos+InPar.BinDisp);
   CutSpat  = median(ImMat(:,MinInd:MaxInd),2);
   
   RegionInd = (floor(LastPeakSpatPos-InPar.PeakSWidth):1:ceil(LastPeakSpatPos+InPar.PeakSWidth)).';
   
   if (~isempty(InPar.Kernel))
      CutSpat  = abs(ifft(fft(CutSpat).*fft(InPar.Kernel,length(CutSpat))));
   end
   
   [PeakSpatPos,PeakH]=find_peak_center(VecSpat,CutSpat,'Method',InPar.PeakMethod,'Ind',RegionInd);
   if (~isnan(PeakSpatPos))
       LastPeakSpatPos = PeakSpatPos;
   end
   Trace.Spat(Ipos)   = PeakSpatPos;
   Trace.Height(Ipos) = PeakH;
end

Trace.HeightErr = sqrt(Trace.Height);

MasterTrace = interp1(InPar.MasterTrace(:,1),InPar.MasterTrace(:,2),Trace.Disp,InPar.InterpMethod');


MasterOffset = (Trace.Spat - nanmedian(Trace.Spat)) - (MasterTrace - nanmedian(MasterTrace));
FlagGoodMT   = (abs(MasterOffset)<InPar.MaxDevMT);
FlagNotEdge  = Trace.Disp>InPar.RemoveEdges & Trace.Disp<(SizeDisp-InPar.RemoveEdges);
FlagGoodErr  = (Trace.HeightErr./Trace.Height)<InPar.MaxRelErr;
FlagGoodDN   = Trace.Height>InPar.MinDN;
FlagGood     = FlagGoodMT & FlagNotEdge & FlagGoodErr & FlagGoodDN;

%[FlagGoodMT,FlagNotEdge,FlagGoodErr]

Trace.Mask   = maskflag_set([],InPar.MaskType,InPar.Bit_BigOffsetMT,~FlagGoodMT,...
                                              InPar.Bit_NearEdge,~FlagNotEdge,...
                                              InPar.Bit_BigErr,~FlagGoodErr);

%--- generate a smooth version of Trace.Spat ---
% interpolate over NaN
Trace.SpatS = medfilt1(interp1_nan(Trace.Disp,Trace.Spat),InPar.BinDisp);

                                          

%-------------------------------------
%--- Fit a polynomial to the trace ---
%-------------------------------------

switch lower(InPar.FitType)
    case 'mt'
        Trace.SpatP = nanmedian(Trace.Spat) + nanmedian(MasterOffset) + (MasterTrace-nanmedian(MasterTrace));
        %plot(Trace.Disp(FlagGood),Trace.Spat(FlagGood),'.')
        %hold on
        %plot(Trace.Disp,Trace.SpatP,'r-')

    case 'poly'
        
        Res = fitgenpoly(Trace.Disp(FlagGood),Trace.Spat(FlagGood),[],...
                         InPar.PolyDeg,InPar.fitgenpolyPar{:});
        Trace.SpatP = polyval(Res.Par,Trace.Disp);
        Poly = Res;
        Trace.Poly  = Poly;
        [Trace.SpatP,~] = polyconf_cov(Poly.Par,Trace.Disp,Poly.Cov);

    otherwise
        error('Unknown FitType option');
end
Trace.Resid    = Trace.Spat - Trace.SpatP;

Flag_BigOffset = abs(Trace.Resid)>InPar.MaxOffset;
Trace.Mask     = maskflag_set(Trace.Mask,InPar.MaskType,...
                              InPar.Bit_BigOffsetFit,Flag_BigOffset);
Trace.RMS      = nanstd(Trace.Resid);
                        

if (InPar.Int),
    % Interactive rejection of points
    %MaxIter    = 1;
    %PlotOption = 1;
    %[XY2,XY1,XY0,Par,Res,IU,INU,GraphicH] = plot_rm(Trace.Disp,Trace.Spat,'b.','polyfit_sc',5,InPar.Deg,InPar.SigClip,MaxIter,PlotOption);

    plot(Trace.Disp(FlagGood),Trace.Spat(FlagGood),'o');
    H = xlabel('Dispersion axis [pix]');
    set(H,'FontSize',16);
    H = ylabel('Spatial axis [pix]');
    set(H,'FontSize',16);
    [Res] = plot_int([],[],...
                     @fitgenpoly,...
                     {[],InPar.PolyDeg,InPar.fitgenpolyPar{:},'Plot','fitonly'},...
                     'FunParInd',InPar.PolyDeg,...
                     'FunBehav','i');

    %waitfor(gcf,'KeyPressFcn','');
    Poly = Res.FunOut{1};
    Trace.Poly  = Poly;
    [Trace.SpatP,~]=polyconf_cov(Poly.Par,Trace.Disp,Poly.Cov);

end



%----------------------------
%--- Extract the spectrum ---
%----------------------------
if (nargout>1)
   ExtractedSpec=spec_extract_2d(ImMat,Trace,varargin{:},'DispDir','x');
end

% delete temp image
if (IsTmp==1)
   delete(Image);
end





