function [SegmentsInfo,AbsCalib]=sedm_wavecalib(varargin)
%--------------------------------------------------------------------------
% sedm_relwavecalib function                                          SEDM
% Description: Find and apply the wavelength calibration to all the
%              spexcells. The function ...
% Input  : * Arbitrary number of pairs of ...,key,val,.. arguments.
%            The following keywords are available:
%            'SI'    - Input SegmentsInfo structure array
%                      as returned by sedm_extract_spexcell.m after it
%                      was run on an ARC image. 
%                      If a string is provided then will assume this is
%                      a mat file containing a structure.
%                      This parameter must be provided.
%            'AbsCalib' - Absolute calibration solution as returned by
%                      sedm_abswavecalib.m. Default is empty matrix.
%                      If empty then attempt to run sedm_abswavecalib.m.
%            'AbsWaveCalibPar' - A cell array of pairs of arguments to pass
%                      to sedm_abswavecalib.m. Default is {} (i.e., see
%                      sedm_abswavecalib.m for default values).
%            'Template' - Use wavelength calibrated template, instead of
%                      solution relative to a specifix spexcell.
%                      This can either a matrix [wave(nm), intensity],
%                      or a mat file name containing [wave(nm), intensity],
%                      or a lamp tempate name {'Hg','HgXe'}.
%                      If empty then use Iref keywork to get ref spectrum.
%                      Default is empty.
%            'Waitbar'- {false, true}. Plot a waitbar showing the fraction
%                      of run time left. Default is false.
%            'SpecField' - The name (string) of the field in the SegmentsInfo
%                       structure array that contains the spectrum in
%                       each spexcell as a function of pixel.
%                       Default is 'SpexSpecFit'.
%            'ScaleVec' - Vector of scales to scan in the first iteration.
%                      Default is (0.8:0.002:1.2).'.
%            'MaxLinesDist' - In the lines matching process, this is the
%                      maximum distance (in nm units) allowed between
%                      matched lines. Default is 6.
%            'Interp1Method' - Interpolation method. See interp1.m for
%                      options. Default is 'linear'.
%            'Interp2Method' - Interpolation method for interpolating
%                      wavelength calibration parameters. See interp2.m for
%                      options. Default is 'cubic'.
%            'MinXCwaveCalib' - Minimum cross-correlation require for
%                      matching. Default is 0.8.
%            'MinNotNaN' - Minimum number of not a NaN in the cross
%                      correlation vector.
%            'List_FluxCorrect' - A cell array of strings of field names
%                      in the SI structure to correct for the d\lambda
%                      effect (e.g., 'SpexSpecFit').
%                      I.e., all the spexcells are divided by
%                      d\lambda. Default is empty cell.
%            'InterpGrid' - Interpolate the grid (true) or use nearest
%                      value (false). Default is false.
% Output : - The input SegmentsInfo structure array with the additional
%            fields related to the wavelength calibration.
%          - The AbsCalib structure (returned by sedm_abswavecalib.m).
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: CArcHgXeSegmentsInfo=sedm_wavecalib('SI',ArcHgXeSegmentsInfo,'AbsCalib',AbsCalib,'waitbar',true)
% Reliable: 2
%--------------------------------------------------------------------------



DefV.SI              = [];
DefV.AbsCalib        = [];
DefV.AbsWaveCalibPar = {};
DefV.Template        = [];  %'HgXe';
DefV.Waitbar         = false;
DefV.SpecField       = 'SpexSpecFit';   % spectrum field to use in SegmentsInfo
DefV.ScaleVec        = (0.8:0.002:1.2).';  %(0.997:0.00005:1.003).';
DefV.MaxLinesDist    = 6;   % [nm]
DefV.Interp1Method   = 'linear';
DefV.Interp2Method   = 'cubic';
DefV.MinXCwaveCalib  = 0.8;
DefV.MinNotNaN       = 5;
DefV.List_FluxCorrect= {};
DefV.InterpGrid      = false;
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});


% read SegmentsInfo
if (isempty(InPar.SI)),
    error('SegmentsInfo (SI) must be provided');
end
if (ischar(InPar.SI)),
    SegmentsInfo = load2(InPar.SI);
else
   SegmentsInfo = InPar.SI;
   InPar = rmfield(InPar,'SI');
end



% AbsCalib
if (isempty(InPar.AbsCalib)),
    AbsCalib = sedm_abswavecalib('SI',SegmentsInfo,InPar.AbsWaveCalibPar{:});
else
    AbsCalib = InPar.AbsCalib;
end

for Igrid=1:1:length(AbsCalib),
    %NN = ~isnan(AbsCalib(Igrid).CalibSpec(:,2));
    [AbsCalib(Igrid).RefLineList]=find_contrast_peaks(AbsCalib(Igrid).CalibSpec(:,1),AbsCalib(Igrid).CalibSpec(:,2),10);
end


Nseg = length(SegmentsInfo);


if (InPar.Waitbar)
   Hwb = waitbar(1,'Fraction of time left: sedm\_wavecalib.m');
end

%
%if (InPar.MultiThred>1),
%   matlabpool(InPar.MultiThred);
%end

MeanX  = [AbsCalib.Xpos].';
MeanY  = [AbsCalib.Ypos].';
BestK  = [AbsCalib.Arc_BestK].';
BestK1 = [AbsCalib.Arc_BestK1].';


for Iseg=1:1:Nseg,
    %Iseg
    %Iseg=1202;
    if (InPar.Waitbar)
       if (floor(Iseg./10)==(Iseg./10)),
          waitbar(1-Iseg./Nseg,Hwb,sprintf('sedm\\_wavecalib.m:  %d/%d done',Iseg,Nseg));
       end
    end
    
    if (InPar.InterpGrid),
        % Interpolate AbsCalib parameters to the spexcell position
        error('InterpGrid=true is not available');
        
        IntBestK = interp2(MeanX,MeanY,BestK,...
                           SegmentsInfo(Iseg).MeanX,SegmentsInfo(Iseg).MeanY,...
                           InPar.Interp2Method,InPar.Interp2Method);
                
        SpecWave = InPar.FunDispLam(BestBestK,SpecPix,BestK1);
        
    else
        % Use nearest AbsCalib 
        % search for nearest calibrated template
        [~,Igrid] = min(Util.Geom.plane_dist(MeanX,MeanY,SegmentsInfo(Iseg).MeanX,SegmentsInfo(Iseg).MeanY));
        % get calibrated spectrum from grid
        CalibTemplateWave = AbsCalib(Igrid).CalibSpec(:,1);
        CalibTemplateInt  = AbsCalib(Igrid).CalibSpec(:,2);
        CalibTemplatePix  = (1:1:length(CalibTemplateWave)).';
    end

    SpecPix = (1:1:length(SegmentsInfo(Iseg).(InPar.SpecField))).';
    SpecInt = SegmentsInfo(Iseg).(InPar.SpecField).';
    %  InPar.ScaleVec = 1.0129;
    Nsv = length(InPar.ScaleVec);
    ResCorr = zeros(Nsv,2);
    for Isv=1:1:Nsv,
        
        Scale = InPar.ScaleVec(Isv);
                     
        %StSpecPix = SpecPix.*Scale;
        %StSpecInt = interp1(SpecPix,SpecInt,StSpecPix,InPar.Interp1Method,0);
        StCalibTemplatePix = CalibTemplatePix.*Scale;
        StCalibTemplateInt = interp1(CalibTemplatePix,CalibTemplateInt,StCalibTemplatePix,InPar.Interp1Method,0);
        SpecInt(isnan(SpecInt))=0; 
        
        
       % plot(CalibTemplatePix,StCalibTemplateInt);
       % hold on
       % plot(SpecInt,'r-')
        
        
        [~,~,Info] = xcorr_fft(StCalibTemplateInt,SpecInt);
        
        ResCorr(Isv,:) = [Info.BestCorr, Info.BestShift];
    end
    
    %plot(InPar.ScaleVec,ResCorr(:,1))
    
    Inn = find(~isnan(ResCorr(:,1)));
    %if (isempty(Inn)),
    if (numel(Inn)<InPar.MinNotNaN),
        % do nothing - bad solution
        SegmentsInfo(Iseg).MaxCorr        = NaN;
        SegmentsInfo(Iseg).Wave_LineResid = NaN;
        SegmentsInfo(Iseg).Wave_LineNum   = 0;
        SegmentsInfo(Iseg).Wave_LineRMS   = NaN;
    else      
        Extram = find_local_extramum(InPar.ScaleVec(Inn),ResCorr(Inn,1));
        Imax = find(Extram(:,3)<0);   % find maxima
        [MaxCorr,MaxInd] = max(Extram(Imax,2));
        BestScale = Extram(Imax(MaxInd),1);
        BestShift = interp1(InPar.ScaleVec,ResCorr(:,2),BestScale,InPar.Interp1Method);
 
        if (isempty(BestScale)),
            BestScale = NaN;
        end
        if (isempty(BestShift)),
            BestShift = NaN;
        end
       % SpecWave = interp1(CalibTemplatePix, CalibTemplateWave, CalibTemplatePix.*BestScale+BestShift,InPar.Interp1Method);
        SpecWave = interp1(CalibTemplatePix, CalibTemplateWave, SpecPix.*BestScale+BestShift,InPar.Interp1Method);
        
%         StCalibTemplatePix = CalibTemplatePix.*BestScale;
%         StCalibTemplateInt = interp1(CalibTemplatePix,CalibTemplateInt,StCalibTemplatePix,InPar.Interp1Method,0); 
%         ShStCalibTemplatePix = StCalibTemplatePix - BestShift.*BestScale;
%         SpecWave = interp1(ShStCalibTemplatePix,CalibTemplateWave,SpecPix,InPar.Interp1Method);

        SegmentsInfo(Iseg).WaveCalib  = SpecWave(:,1).';
        SegmentsInfo(Iseg).MaxCorr    = MaxCorr;
        SegmentsInfo(Iseg).BestScale  = BestScale;
        SegmentsInfo(Iseg).BestShift  = BestShift;
        SegmentsInfo(Iseg).Igrid      = Igrid;
       
        % coorect flux for d\lambda effect
        if (~isempty(InPar.List_FluxCorrect)),
            DLambda = zeros(size(SegmentsInfo(Iseg).WaveCalib));
            DLambda = diff(SegmentsInfo(Iseg).WaveCalib);
            DLambda = [DLambda(1); DLambda];
            for Ifluxc=1:1:numel(InPar.List_FluxCorrect),
                SegmentsInfo(Iseg).(InPar.List_FluxCorrect{Ifluxc})./DLambda;
            end
        end
        
        
        % estimate solution RMS using lines
        %plot(SpecWave,SpecInt)
       
        [LineList]=find_contrast_peaks(SpecWave,SpecInt,10);
        RefLineList = AbsCalib(Igrid).RefLineList;

        LinesCounter = 0;
        %clear MatchedLines;
        MatchedLines = zeros(size(RefLineList,1),2).*NaN;                                            

        for Irll=1:1:size(RefLineList,1),
            Dist = RefLineList(Irll,1) - LineList(:,1);
            Ill = find(abs(Dist)<InPar.MaxLinesDist);
            if (~isempty(Ill)),
                LinesCounter = LinesCounter + 1;    
                if (length(Ill)>1),
                    % if more than one line - choose highest (rather
                    % than nearest).
                    [~,MaxIll] = max(LineList(Ill,2));
                    Ill = Ill(MaxIll);
                end
                MatchedLines(Irll,:) = [RefLineList(Irll,1), LineList(Ill,1)];
            end
        end
        %MatchedLines
        Ilnn = find(~isnan(MatchedLines(:,2)));
        SegmentsInfo(Iseg).Wave_LineResid = MatchedLines(Ilnn,1) - MatchedLines(Ilnn,2);
        SegmentsInfo(Iseg).Wave_LineNum   = length(Ilnn);
        SegmentsInfo(Iseg).Wave_LineRMS   = std(SegmentsInfo(Iseg).Wave_LineResid);
        
    end
    %MaxCorr
    %BestShift
    %BestScale
    %plot(SpecWave,SpecInt);
    %hold on;
    %plot(CalibTemplateWave,CalibTemplateInt,'r-');
    %drawnow
        
end

if (InPar.Waitbar)
   close(Hwb);
end
       
    
