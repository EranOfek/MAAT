function SegmentsInfo=sedm_extract_spexcell(varargin)
%--------------------------------------------------------------------------
% sedm_extract_spexcell function                                      SEDM
% Description: Given a SEDM SegmentsInfo structure returned by
%              sedm_generate_spexcell_segmentation.m and a science or
%              calibration image, search for the exact position and slope
%              of the trace in each spexcell using the SegmentsInfo as a
%              first guess.
% Input  : * Arbitrary number of pairs of ...,key,val,.. arguments.
%            The following keywords are available:
%            'SI'    - Input SegmentsInfo structure as returned by
%                      sedm_refine_slopepos.m.
%                      If a string is provided then will assume this is
%                      a mat file containing a structure named SegmentsInfo.
%                      This parameter must be provided.
%            'ScienceImage' - The science/calibration image from which to
%                      extract the spectrum. This can be a matrix, string
%                      cell, structure or any of the inputs allowed by
%                      read2sim.m
%                      This parameter must be provided.
%            'Gain'  - The CCD Gain (scalar) or a the CCD Gain header
%                      keyword (string). Default is 'GAIN'.
%                      All the output spexcells will returned in e- units.
%                      If empty, don't apply Gain.
%            'OffsetType' - {'measured','smooth'} indicate which trace
%                      offset estimate to use. The one directly measured
%                      in each spexcell (i.e., SegmentsInfo.MeasuredOffset)
%                      or the fitted polynomial surface
%                      (i.e., SegmentsInfo.PredictedOffset).
%                      Default is 'smooth'.
%            'SemiWidthExtractionBlock' - Semi width of the extracted
%                      trace around the spexcell along the position
%                      dimension. The total width of the extraction block
%                      will be this parameter*2 +1. Default is 4.
%            'SemiWidthMeasuredBlock' - Semi width of the measured
%                      extracted trace around the spexcell along the
%                      position dimension. The total width of the
%                      extraction block will be this parameter*2 +1.
%                      Default is 3.
%            'StepExtractionBlock' - Interpolation step size of the
%                      extracted trace around the spexcell along the
%                      position dimension. Default is 1. If values smaller
%                      than 1 are used then the spectra will be over
%                      sampled.
%            'XstartOffset' - slope start point relative to segment start.
%                      Default is 30.
%            'XendOffset' - slope end point relative to segment end.
%                      Default is 130.
%            'Interp2Method' - 2D interpolation method (see interp2.m
%                      for option). Default is 'cubic'.
%            'SemiWidthProf' - Semi width of the section from which the
%                      spexcell profile will be extracted around the
%                      maximum intensity. Default is 15 pix.
%            'SubSpexcellBack' - Subtract background from each spexcell.
%                      Options are:
%                      'none' - Do not subtract background
%                      'all'  - Subtract background where the background
%                               is estimated from the 3-sigma lower
%                               percentile of the entire image.
%                      'spex' - Subtract background where the background
%                               is estimated from the edges of each
%                               spexcell block.
%                      'spexfit' - Fit a line for each dispersion position
%                               in each spexcell.
%            'SpexFitMF' -Size of median filtering for the 'spexfit'
%                       option. Default is 20. If 0 don't median filter.
%            'IgnoreEdge' - Number of pixels in edge of block to ignore in
%                      the spectrum fitting. Default is 2 (i.e., with
%                      default parameters will use pixels 3:7).
%            'MultiThred' - Number of processors threds to use.
%                      Default is 1.
%                      Run time 3.2s. NOT AVAILABLE YET.
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
%            .Back           - The background in ADU subtracted from each
%                              spexcell.
%          - A structure with the parameters of the 2-D best fit surface
%            fitted to the offset as a function of position.
%            .Res    - the results from fit_2d_polysurface.m
%            .FunStr - A string containing the 2D fitted function.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% SegmentsInfo=sedm_extract_spexcell('SI',SegmentsInfo,'ScienceImage',List.LampHg{1});
% SegmentsInfo=sedm_extract_spexcell('SI',SegmentsInfo,'ScienceImage','b_ifu20130809_20_44_38.fits','MultiThred',4);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.fit.*

ImageField    = 'Im';   % read2sim structure definitions

DefV.SI                        = [];
DefV.ScienceImage              = [];
DefV.Gain                      = 'GAIN';
DefV.OffsetType                = 'smooth';   % {'measured','smooth'}
DefV.SemiWidthExtractionBlock  = 4;
DefV.SemiWidthMeasuredBlock    = 3;
DefV.StepExtractionBlock       = 1;
DefV.XstartOffset              = 30;   
DefV.XendOffset                = 130;
DefV.Interp2Method             = 'cubic';
DefV.SemiWidthProf             = 15;
DefV.SubSpexcellBack           = 'all';   %{'spex','all','non'}
DefV.SpexFitMF                 = 20;
DefV.IgnoreEdge                = 2;
DefV.MultiThred                = 1;
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

if (isempty(InPar.SI)),
    error('SegmentsInfo structure (SI) must be provided');
end
if (isempty(InPar.ScienceImage)),
    error('ScienceImage must be provided');
end

if (ischar(InPar.SI)),
    SegmentsInfo = load2(InPar.SI);
else
   SegmentsInfo = InPar.SI;
   InPar = rmfield(InPar,'SI');
end

SIM = read2sim(InPar.ScienceImage);
ScienceImage = SIM.(ImageField);
%InPar = rmfield(InPar,'ScienceImage');

if (~isempty(InPar.Gain)),
    if (ischar(InPar.Gain)),
        % get Gain from header keyword
        if (ischar(InPar.ScienceImage)),
            [KeywordVal]=get_fits_keyword(InPar.ScienceImage,{InPar.Gain});
            InPar.Gain = KeywordVal{1};
        else
            % can't read Gain - set it to empty
            warning('Cannot read Gain from header - set Gain to empty');
            InPar.Gain = [];
        end
    end
end
            
           
% image size vectors
VecX        = (1:1:size(ScienceImage,2)).';
VecY        = (1:1:size(ScienceImage,1)).';


% Construct the appropriate trace parameters for each spexcell
AllPar1 = [SegmentsInfo.Par1]';  % get fitted traces in segmentation image
switch lower(InPar.OffsetType)
    case 'measured'
        SpexOffset = [SegmentsInfo.MeasuredOffset].';
        % in rare cases MeasuredOffset is NaN
        % in this case replace it with PredictedOffset
        InSO = find(isnan(SpexOffset));
        SpexOffset(InSO) = [SegmentsInfo(InSO).PredictedOffset].';
        
    case 'smooth'
        SpexOffset = [SegmentsInfo.PredictedOffset].';
    otherwise
        error('Unknown OffsetType option');
end
AllPar1(:,2) = AllPar1(:,2) + SpexOffset;



% some statistics for background subtraction
ErrCL = err_cl(ScienceImage(:),0.9973);


%--------------------------------------------------------
%--- Remap and extract spexcell spectrum around trace ---
%--------------------------------------------------------
BlockRelY = (-InPar.SemiWidthExtractionBlock:InPar.StepExtractionBlock:InPar.SemiWidthExtractionBlock);
Nby       = length(BlockRelY);
Nseg      = length(SegmentsInfo);



if (InPar.MultiThred>1),
   matlabpool(InPar.MultiThred) 
end
%parfor Iseg=1:1:Nseg,
for Iseg=1:1:Nseg,
    Step   = 1;
    XI      = (SegmentsInfo(Iseg).MinX-InPar.XstartOffset:Step:SegmentsInfo(Iseg).MaxX+InPar.XendOffset).'-SegmentsInfo(Iseg).MeanX;
   
    YI      = polyval(AllPar1(Iseg,:).',XI);
    XI      = XI + SegmentsInfo(Iseg).MeanX;
    
    % X and Y position of the extrcated block center
    SegmentsInfo(Iseg).BlockCenterX = XI(InPar.SemiWidthExtractionBlock+1,:);
    SegmentsInfo(Iseg).BlockCenterY = YI(InPar.SemiWidthExtractionBlock+1,:);
    
    Slope             = AllPar1(Iseg,1);
    CosArcTanSlope    = sqrt(1 + Slope.^2);
    CosArcTanInvSlope = sqrt(1 + Slope.^(-2));
    %GridX             = zeros(Nby,length(XI));
    %GridY             = GridX;
   
    % calculate XI and YI which are parallel to the original XI YI and
    % are shifted upward/downward by an given offset
    Delta = BlockRelY ./ CosArcTanSlope;
    GridX = bsxfun(@minus,XI,sign(Slope).*Delta./CosArcTanInvSlope).';
    GridY = bsxfun(@plus,YI,(bsxfun(@minus,XI,GridX.'))./Slope).';     
    
    % interpolate spexcell to new uniform and parallel grid    
    
    SpexSpecBlock = interp2fast(VecX,VecY,ScienceImage,...
                                GridX,GridY,InPar.Interp2Method);

    %if (~isempty(InPar.Gain)),
    %    SpexSpecBlock = SpexSpecBlock.*InPar.Gain;
    %    SegmentsInfo(Iseg).Units = 'e-';
    %else
    %    SegmentsInfo(Iseg).Units = 'ADU';
    %end
    
    
    SegmentsInfo(Iseg).SpexSpecBlock = SpexSpecBlock;
    
    % Extract 1d spectrum for each spex
    % currently use a simplified approach
    % This section need to be developed!
    SegmentsInfo(Iseg).SpexSpecCenter = nansum(SegmentsInfo(Iseg).SpexSpecBlock(4:6,:),1);
    [~,MaxInd] = max(SegmentsInfo(Iseg).SpexSpecCenter);
    
    % Numerical profile of spexcell
    SegmentsInfo(Iseg).Profile = nanmedian(SegmentsInfo(Iseg).SpexSpecBlock,2);
    % central profile around max
    CInd   = Util.array.index_outofbound((MaxInd-InPar.SemiWidthProf:1:MaxInd+InPar.SemiWidthProf),size(SegmentsInfo(Iseg).SpexSpecBlock,2));
    SegmentsInfo(Iseg).ProfileMax = nanmedian(SegmentsInfo(Iseg).SpexSpecBlock(:,CInd),2);
   
    %--- subtract background from spexcell ---
    % background subtracted SpexSpecBloc
    
    switch lower(InPar.SubSpexcellBack)
        case 'none'
            % no background subtraction in spex
            Back = 0;
            SegmentsInfo(Iseg).BS_SpexSpecBlock = SegmentsInfo(Iseg).SpexSpecBlock - Back;
            SegmentsInfo(Iseg).BS_Profile       = SegmentsInfo(Iseg).Profile       - Back;
            SegmentsInfo(Iseg).BS_ProfileMax    = SegmentsInfo(Iseg).ProfileMax    - Back;
        case 'all'
            % Background is estimated from the global properties
            % of the image using the lower 3sigma percentile
            Back = ErrCL(1,1);
            SegmentsInfo(Iseg).BS_SpexSpecBlock = SegmentsInfo(Iseg).SpexSpecBlock - Back;
            SegmentsInfo(Iseg).BS_Profile       = SegmentsInfo(Iseg).Profile       - Back;
            SegmentsInfo(Iseg).BS_ProfileMax    = SegmentsInfo(Iseg).ProfileMax    - Back;
        case 'spex'
            % Background is estimated from the edges of the spexcell
            Back = median([SegmentsInfo(Iseg).SpexSpecBlock(1,:), SegmentsInfo(Iseg).SpexSpecBlock(Nby,:)]);
            SegmentsInfo(Iseg).BS_SpexSpecBlock = SegmentsInfo(Iseg).SpexSpecBlock - Back;
            SegmentsInfo(Iseg).BS_Profile       = SegmentsInfo(Iseg).Profile       - Back;
            SegmentsInfo(Iseg).BS_ProfileMax    = SegmentsInfo(Iseg).ProfileMax    - Back;
        case 'spexfit'
            % Estimate background by fitting a line in each spexcell
            % disperssion position. The line is fitted only useing the
            % edges.
            IndSpex = [1;2; Nby-1; Nby];
            Hspex   = [ones(4,1), IndSpex];
            Yspex   = SegmentsInfo(Iseg).SpexSpecBlock(IndSpex,:);
            ParSpex = Hspex\Yspex;
            
            HspexA  = [ones(Nby,1), (1:1:Nby).'];
            Back    = HspexA*ParSpex;
             % median filter Back
            if (InPar.SpexFitMF>0),
                Back = medfilt1(Back,InPar.SpexFitMF,[],2);
            end
            SegmentsInfo(Iseg).BS_SpexSpecBlock = SegmentsInfo(Iseg).SpexSpecBlock - Back;
            SegmentsInfo(Iseg).BS_Profile       = nanmedian(SegmentsInfo(Iseg).SpexSpecBlock - Back,2);
            SegmentsInfo(Iseg).BS_ProfileMax    = nanmedian(SegmentsInfo(Iseg).SpexSpecBlock(:,CInd) - Back(:,CInd),2);
            
        otherwise
            error('Unknown SubSpexcellBack option');
    end
    SegmentsInfo(Iseg).Back = Back;
    SegmentsInfo(Iseg).MedBack = median(Back(:));   % median background in spexcell
    
    % profile normalization    
    SegmentsInfo(Iseg).ProfileNorm       = trapz(SegmentsInfo(Iseg).BS_ProfileMax);
    % background subtracted normalized profile
    SegmentsInfo(Iseg).NormalizedProfile = SegmentsInfo(Iseg).BS_ProfileMax./SegmentsInfo(Iseg).ProfileNorm;  
    
    % simultanous least square fitting of the entire spectrum
    Y = SegmentsInfo(Iseg).BS_SpexSpecBlock(1+InPar.IgnoreEdge:1:Nby-InPar.IgnoreEdge,:);
    H = SegmentsInfo(Iseg).NormalizedProfile(1+InPar.IgnoreEdge:1:Nby-InPar.IgnoreEdge);
    Par = inv(H.'*H)*H.'*Y;
    %Par = inv(SegmentsInfo(Iseg).NormalizedProfile.'*SegmentsInfo(Iseg).NormalizedProfile)*SegmentsInfo(Iseg).NormalizedProfile.'*SegmentsInfo(Iseg).BS_SpexSpecBlock;
    SegmentsInfo(Iseg).SpexSpecFit        = Par;
    % quick error estimate
    % median(Back,1) is needed for the spexfot option
    SegmentsInfo(Iseg).SpexSpecFitErr     = sqrt(SegmentsInfo(Iseg).SpexSpecFit + median(Back,1)).*sqrt(InPar.Gain);
    SegmentsInfo(Iseg).SpexSpecFitResid   = Y - H*Par;
    SegmentsInfo(Iseg).SpexSpecFitChi2    = sum(SegmentsInfo(Iseg).SpexSpecFitResid.^2./SegmentsInfo(Iseg).SpexSpecBlock(1+InPar.IgnoreEdge:1:Nby-InPar.IgnoreEdge,:),1);
    
    
    
    if (isempty(find(isnan(H),1))),
        SegmentsInfo(Iseg).FitFlag = cond(H.'*H);
    else
        SegmentsInfo(Iseg).FitFlag = 0;
    end
    
    %for Ipix=1:1:length(SegmentsInfo(Iseg).SpexSpec),
    %    % fit the profile
    %    if (isnan(SegmentsInfo(Iseg).ProfileNorm)),
    %        % need to deal with NaNs
    %    else
    %       [Par,ParErr]=lscov(SegmentsInfo(Iseg).NormalizedProfile,...
    %                          SegmentsInfo(Iseg).SpexSpecBlock(:,Ipix),...
    %                          SegmentsInfo(Iseg).Profile);
    %       SegmentsInfo(Iseg).SpexSpec(Ipix)        = Par;
    %       SegmentsInfo(Iseg).SpexSpecErr(Ipix)     = ParErr;
    %       SegmentsInfo(Iseg).SpexSpecResid(:,Ipix) = SegmentsInfo(Iseg).SpexSpecBlock(:,Ipix) - ...
    %                                                  SegmentsInfo(Iseg).NormalizedProfile*Par;
    %       SegmentsInfo(Iseg).SpexSpecChi2(Ipix)    = sum(SegmentsInfo(Iseg).SpexSpecResid(:,Ipix).^2./SegmentsInfo(Iseg).SpexSpecBlock(:,Ipix));
    %    end
    %end
    
    
end
if (InPar.MultiThred>1),
    matlabpool('close');
end






