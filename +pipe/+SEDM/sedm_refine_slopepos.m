function [SegmentsInfo,OffsetSurface]=sedm_refine_slopepos(varargin)
%--------------------------------------------------------------------------
% sedm_refine_slopepos function                                       SEDM
% Description: Given a SEDM SegmentsInfo structure returned by
%              sedm_generate_spexcell_segmentation.m and a science or
%              calibration image, search for the exact position and slope
%              of the trace in eac spexcell using the SegmentsInfo as a
%              first guess.
% Input  : * Arbitrary number of pairs of ...,key,val,.. arguments.
%            The following keywords are available:
%            'Input' - A single bias subtracted science or calibration
%                      image in which to refine the slopes in each
%                      spexcell. This can be:
%                      a string containing image name,
%                      or a matrix containing the image, or a structure
%                      in which the .Im fields contain the image matrix.
%                      There is no default value for this parameter,
%                      and it must be provided.
%            'SI'    - SegmentsInfo structure as returned by
%                      sedm_generate_spexcell_segmentation.m.
%                      If a string is provided then will assume this is
%                      a mat file containing a structure named SegmentsInfo.
%                      Default is empty matrix.
%                      If empty matrix will attempt to generate
%                      SegmentsInfo by calling
%                      sedm_generate_spexcell_segmentation.m with default
%                      parameters.
%            'SlopeType' - {'measured','smooth'} indicate which trace slope
%                      estimate to use. The one directly measured in each
%                      spexcell (i.e., SegmentsInfo.Par1(1)) or
%                      the fitted polynomial surface
%                      (i.e., SegmentsInfo.PredictedSlope).
%                      Default is 'smooth'.
%            'OffsetYvec' - Y offset to scan while looking for the best
%                      trace position (by maximizing the integral along
%                      the trace).
%            'Interp2Method' - 2D interpolation method (see interp2.m
%                      for option). Default is 'cubic'.
%            'XstartOffset' - slope start point relative to segment start.
%                      Default is 30.
%            'XendOffset' - slope end point relative to segment end.
%                      Default is 130.
%            'PolyXslope' - Orders of polynomial X axis orders to use
%                      in the fit of the 2-D surface to the slope as
%                      a function of position.
%                      Default is [1 2].
%            'PolyYslope' - Orders of polynomial Y axis orders to use
%                      in the fit of the 2-D surface to the slope as
%                      a function of position.
%                      Default is [1 2].
%            'PolyXYslope' - Orders of polynomial XY axis orders to use
%                      in the fit of the 2-D surface to the slope as
%                      a function of position.
%                      Default is [1 1].
%            'MultiThred' - Number of processors threds to use.
%                      Default is 1.
%                      Run time is 36, 41 and 50 s for 8, 4 and 1 threds,
%                      respectively.
% Output : - Structure array of segmentation and slope info.
%            Each elemnt contains a structure of information about
%            a specific segment.
%            The structure includes all the fields returned by
%            sedm_generate_spexcell_segmentation. In addition the
%            following fields are available:
%            .MeasuredOffset - The measured Y offset between the
%                              slope available from the segmentation 
%                              image and the slope measured in the
%                              current image.
%            .PredictedOffset - Predicted offset by fitting a polynomial
%                              surface to the MeasuredOffset offset
%                              parameter, as a function of x/y position.
%          - A structure with the parameters of the 2-D best fit surface
%            fitted to the offset as a function of position.
%            .Res    - the results from fit_2d_polysurface.m
%            .FunStr - A string containing the 2D fitted function.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [SegmentsInfo,OffsetSurface]=sedm_refine_slopepos('Input','b_ifu20130809_05_26_52.fits','SI',SegmentsInfo);
%          [SegmentsInfo,OffsetSurface]=sedm_refine_slopepos('Input','b_ifu20130809_05_26_52.fits');
%          [SegmentsInfo,OffsetSurface]=sedm_refine_slopepos('Input','b_ifu20130809_20_44_38.fits','SI',SegmentsInfo);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.fit.*

DefV.Input                  = [];
DefV.SI                     = [];
DefV.SlopeType              = 'smooth';   % {'measured','smooth'}
DefV.OffsetYvec             = (-5:0.5:5).';
DefV.Interp2Method          = 'cubic';
DefV.XstartOffset           = 0;  %30;   
DefV.XendOffset             = 0;  %130;
DefV.PolyXoffset            = [1 2];
DefV.PolyYoffset            = [1 2];
DefV.PolyXYoffset           = [1 1];
DefV.MultiThred             = 1;
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

if (isempty(InPar.Input)),
    error('Input keyword argument must be provided');
end

if (isempty(InPar.SI)),
    [SegmentsInfo]=sedm_generate_spexcell_segmentation;
else
    if (ischar(InPar.SI)),
        load(InPar.SI);
        if (exist('SegmentsInfo','var')~=1),
            error('Can not find SegmentsInfo structure in %s',InPar.SI);
        end
    else
        % InPar.SI is a structure.
        SegmentsInfo = InPar.SI;
    end
end
Nseg = length(SegmentsInfo);


if (ischar(InPar.Input)),
    ScienceImage     = fitsread(InPar.Input);
elseif (isstruct(InPar.Input)),
    ScienceImage = InPar.Input.Im;
elseif (isnumeric(InPar.Input)),
    ScienceImage = InPar.Input;
else
    error('Unknown Input data type');
end


% define which slope to use
AllPar1 = [SegmentsInfo.Par1].';
switch lower(InPar.SlopeType)
    case 'measured'
        % do nothing - use measured slope and offsets
    case 'smooth'
        AllPar1(:,1) = [SegmentsInfo(1:1:Nseg).PredictedSlope].';
    otherwise
        error('Unknown SlopeType option');
end

% image size vectors
VecX        = (1:1:size(ScienceImage,2)).';
VecY        = (1:1:size(ScienceImage,1)).';

%--------------------------------------------------------------------
%--- Find best Y-axis offset of trace relative to reference trace ---
%--------------------------------------------------------------------
% in each spexcell estimate the Y-axis offset between spectrum position and
% first guess spectrum position

if (InPar.MultiThred>1),
   matlabpool(InPar.MultiThred) 
end
for Iseg=1:1:Nseg,
%parfor Iseg=1:1:Nseg,
    %SegmentsInfo(Iseg)
    %Iseg
    Step   = 1;
    XI      = (SegmentsInfo(Iseg).MinX-InPar.XstartOffset:Step:SegmentsInfo(Iseg).MaxX+InPar.XendOffset).'-SegmentsInfo(Iseg).MeanX;
    YI      = polyval(AllPar1(Iseg,:).',XI);
    XI      = XI + SegmentsInfo(Iseg).MeanX;
    % find best shift in Y axis between segmentation trace and spectrum
    % trace for specific spexcell
    
    % select maxima as a best offset solution
    % check if maxima is not too far away...
    SumLight   = zeros(length(InPar.OffsetYvec),1);
    for Ioffset=1:1:length(InPar.OffsetYvec),  
       InterpTrace = interp2fast(VecX,VecY,ScienceImage,XI,YI+InPar.OffsetYvec(Ioffset),InPar.Interp2Method);
       SumLight(Ioffset) = nansum(InterpTrace);  
    end
    
    Extram = find_local_extramum(InPar.OffsetYvec,SumLight);
    if (isempty(Extram)),
        % no local solution was found
        MeasuredOffset = NaN;
    else
        Extram = Extram(Extram(:,3)<0,:);
        if (isempty(Extram)),
            MeasuredOffset = NaN;
        else
           [~,MaxI] = max(Extram(:,2));       
           MeasuredOffset = Extram(MaxI,1);
        end
    end
    
    %[~,MaxI]   = max(SumLightI);
    SegmentsInfo(Iseg).MeasuredOffset = MeasuredOffset;
    
    %plot(InPar.OffsetYvec,SumLight)
    %plot(InPar.OffsetYvecI,SumLightI)
    
end
if (InPar.MultiThred>1),
   matlabpool('close');
end



% fit a surface to MeasuredOffset parameter
Xa = [SegmentsInfo.MeanX].';
Ya = [SegmentsInfo.MeanY].';
MeasuredOffset = [SegmentsInfo.MeasuredOffset].';

[Res,FunStr]=fit_2d_polysurface(Xa,Ya,MeasuredOffset,[],'X',InPar.PolyXoffset,...
                                                        'Y',InPar.PolyYoffset,...
                                                        'XY',InPar.PolyXYoffset);
OffsetSurface.Res    = Res;
OffsetSurface.FunStr = FunStr;
     
%ParOffsetSurface = Res.Par;
PredictedOffset  = Res.PredZ;
%Resid            = Res.Resid;
PredOffCell = num2cell(PredictedOffset);
[SegmentsInfo(1:1:Nseg).PredictedOffset] = deal(PredOffCell{:});



%InPar.OffsetType               = 'smooth';       % {'measured','smooth'}
%
%switch lower(InPar.OffsetType)
%    case 'measured'
%        SpexOffset = [SegmentsInfo.MeasuredOffset].';
%    case 'smooth'
%        SpexOffset = [SegmentsInfo.PredictedOffset].';
%    otherwise
%        error('Unknown OffsetType option');
%end
%AllPar1(:,2) = AllPar1(:,2) + SpexOffset;




