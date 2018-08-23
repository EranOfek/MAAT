function [SegmentsInfo,SegMap,SlopeSurface]=sedm_generate_spexcell_segmentation(varargin)
%--------------------------------------------------------------------------
% sedm_generate_spexcell_segmentation function                        SEDM
% Description: Use SEDM dome flat images to generatre segmentation map
%              of the individual spexcells (i.e., each segment is marked
%              by a serial number).
% Input  : * Arbitrary number of pairs of ...,key,val,.. arguments.
%            The following keywords are available:
%            'Input' - Bias subtracted dome-flat input images.
%                      String containing images template, or file
%                      containing list of images or cell array containing
%                      image names (see create_list.m for details).
%                      The list can contain either the dome flat
%                      images or all the images. In the later case, the
%                      user should supply information on how to identify
%                      the flat images (see below).
%                      Default is 'b_ifu*.fits'. If empty use default.
%            'IsListFF'- {'y'|'n'}, indicating if the 'Input' list is
%                      a list of only dome flat ('y') or all images ('n').
%                      If this is a list of all images then the program
%                      will attempt to look for the dome flats among
%                      the images. Default is 'y'.
%            'HeadKeys'- Cell array of header keywords which are required
%                      in order to look for dome flat images.
%                      Default is {'OBJECT','NAME'}.
%            'HeadKeysVal' - Cell array of the substrings values which each
%                      one of the keywords listed in 'HeadKeys' can take
%                      place in case the image is dome flat.
%                      Default is {'dome','STOW:Flats'}.
%            'SaturationLevel' - Saturation level of images.
%                      This can be either a scalar or a string of an
%                      image header keyword indicating the saturation
%                      level. Default is 55000.
%                      This parameter will be used to select flat images
%                      only if IsListFF is 'n'.
%            'IFU_FlatImageName' - Name of output combined flat image name.
%                      Default is 'Flat.fits'.
%            'IFU_FlatStDImageName' - Name of output flat StD image.
%                      Default is 'FlatStD.fits'.
%            'ThresholdVal' - Thresholding for segmentation.
%                      Default is 2.7.
%            'SegmentationImage' - Output segmentation image name.
%                      Default is 'SegmentationMap.fits'.
%            'FitSlope' - {'y' | 'n'} fit slopes to each segment.
%                      Default is 'y'.
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
%                      Default is 1. NOT AVAILABLE.
% Output : - Structure array of segmentation info.
%            Each elemnt contains a structure of information about
%            a specific segment. The following fields are available:
%            .SegArea         - segment area [pix].
%            .MeanX           - Mean X position of segments pixels.
%            .MeanY           - Mean Y position of segments pixels.
%            .MinX            - Minimum X position of segments pixels.
%            .MaxX            - Maximum X position of segments pixels.
%            .MinY            - Minimum Y position of segments pixels.
%            .MaxY            - Maximum Y position of segments pixels.
%            .RangeX          - Range of segments X position.
%            .RangeY          - Range of segments Y position.
%            .FlagGood        - Flag indicating if segment best fit
%                               first and second order polynomials are
%                               deviating from each other.
%            .Par1            - Parameters of 1st order polynomial fit
%                               to segment [slope, intersection].
%            .ParErr1         - Errors in Par1.
%            .Par2            - Parameters of 2nd order polynomial fit
%                               to segment [Slope2, slope, intersection].
%            .ParErr2         - Errors in Par2.
%            .PredictedSlope  - Slope predicted by fitting a 2-D surface
%                               to the slope as a function of position.
%            .BadSegmentFlag  - A flag indicating if segment length is
%                               longer than 170 pixels.
%          - Matrix of segmentation map.
%          - A structure with the parameters of the 2-D best fit surface
%            fitted to the slope as a function of position.
%            .Res    - the results from fit_2d_polysurface.m
%            .FunStr - A string containing the 2D fitted function.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [SegmentsInfo,SegMap]=sedm_generate_spexcell_segmentation('Input',List.DomeFlat);
%          [SegmentsInfo,SegMap]=sedm_generate_spexcell_segmentation('IsListFF','n');
%          [SegmentsInfo,SegMap,SlopeSurface]=sedm_generate_spexcell_segmentation;
% Reliable: 2
%--------------------------------------------------------------------------
import Util.fit.*

DefV.Input                  = 'b_ifu*.fits';
DefV.IsListFF               = 'y';
DefV.HeadKeys               = {'OBJECT','NAME'};
DefV.HeadKeysVal            = {'dome','STOW:Flats'};
DefV.SaturationLevel        = 55000;
DefV.IFU_FlatImageName      = 'Flat.fits';
DefV.IFU_FlatStDImageName   = 'FlatStD.fits';
% segmentation parameters
DefV.ThresholdVal           = 2.7;  % good for DomeFlat
DefV.SegmentationImage      = 'SegmentationMap.fits';   % output segmentation image name
% slope fitting
DefV.FitSlope               = 'y';  % {'y' | 'n'}
DefV.MultiThred             = 1;   % number of threds to use in slope fitting
DefV.XstartOffset           = 30;   
DefV.XendOffset             = 130;
DefV.PolyXslope             = [1 2];
DefV.PolyYslope             = [1 2];
DefV.PolyXYslope            = [1 1];

%DefV.SlopeType = 'smooth';   % {'measured','smooth'}
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

[~,ListCell] = create_list(InPar.Input,NaN);

switch lower(InPar.IsListFF)
    case 'y'
        % do nothing
        List.DomeFlat = ListCell;
    case 'n'
        % get some info on images.
        if (ischar(InPar.SaturationLevel)),
            InPar.HeadKeys{end+1} = InPar.SaturationLevel;
            SL_Header = true;
        else
            SL_Header = false;
        end 
            
        %[~,KeywordS] = mget_fits_keyword(ListCell,InPar.HeadKeys);
        [~,~,KeywordS] = fits_mget_keys(ListCell,InPar.HeadKeys);
        Stat         = imstat_fits(ListCell);

        % identify dome flats
        if SL_Header
            IFU_SaturationLevel = [KeywordS.(InPar.SaturationLevel)];
        else
            IFU_SaturationLevel = InPar.SaturationLevel.*ones(1,length(Stat));
        end
            
        % select dome flats
        % cellfun(@num2str,{KeywordS.NAME},'UniformOutput',false) is
        % required in order to convert numbers to strings
        FlagDomeFF = (Util.cell.isempty_cell(strfind(cellfun(@num2str,{KeywordS.NAME},'UniformOutput',false),InPar.HeadKeysVal{2}))==0) & ...
                     (Util.cell.isempty_cell(strfind(cellfun(@num2str,{KeywordS.OBJECT},'UniformOutput',false),InPar.HeadKeysVal{1}))==0) & ...
                     [Stat.Max]<IFU_SaturationLevel;

        List.DomeFlat = ListCell(FlagDomeFF);
    otherwise
        error('Unknown IsListFF option');
end
 
%---------------------
%--- Flat fielding ---
%---------------------
% generate super dome flat image
flatset_fits([],List.DomeFlat,[],InPar.IFU_FlatImageName,...
             'FlatCor','n',...
             'StD',InPar.IFU_FlatStDImageName);

% read flat image
FlatMat = fitsread(InPar.IFU_FlatImageName);

%-------------------------------------
%--- Generate the segmentation map ---
%-------------------------------------
% thresholding and generating spexcell mask image                                                    
OutMat=imreplace_fits('temp_seg1.fits',InPar.IFU_FlatImageName,'Range',[-Inf InPar.ThresholdVal;InPar.ThresholdVal Inf],'Value',[0;1]);

% smoothing to avoid disconnected areas
ThreshMap = conv2(OutMat{1},ones(3,25),'same');
fitswrite(ThreshMap,'temp_seg2.fits');

% thresholding again
imreplace_fits('temp_seg3','temp_seg2.fits','Range',[-Inf 1;1 Inf],'Value',[0;1]);
SegmentationMat = fitsread('temp_seg3');

% watershed transform
WaterShed = watershed(-SegmentationMat);
%fitswrite(WaterShed,'ttt.fits');
SegMap = uint16(WaterShed).*uint16(SegmentationMat);
fitswrite(SegMap,InPar.SegmentationImage);



%-----------------------------------------------------------------
%--- Fit a slope to the pixels (in flat image) of each segment ---
%-----------------------------------------------------------------
switch lower(InPar.FitSlope)
    case 'n'
         SegmentsInfo = [];
         SlopeSurface = [];
    case 'y'
        % go over all segments and extract spectra
        Nseg = max(max(SegMap));
        
        
        VecX        = (1:1:size(SegMap,2)).';
        VecY        = (1:1:size(SegMap,1)).';
        [MatX,MatY] = meshgrid(VecX,VecY);

        %tic;
        %if (InPar.MultiThred>1),
        %   matlabpool(InPar.MultiThred) 
        %end
        SegmentsInfo = struct('SegArea',cell(Nseg,1),...
                              'MeanX',cell(Nseg,1),...
                              'MeanY',cell(Nseg,1),...
                              'MinX',cell(Nseg,1),...
                              'MaxX',cell(Nseg,1),...
                              'MinY',cell(Nseg,1),...
                              'MaxY',cell(Nseg,1),...
                              'RangeX',cell(Nseg,1),...
                              'RangeY',cell(Nseg,1),...
                              'FlagGood',cell(Nseg,1),...
                              'Par1',cell(Nseg,1),...
                              'Par2',cell(Nseg,1),...
                              'PredictedSlope',cell(Nseg,1),...
                              'BadSegmentFlag',cell(Nseg,1));
        
        %parfor Iseg=1:1:Nseg,
        for Iseg=1:1:Nseg,   % use this line if yo don't have the parallel toolbox - use
            IndPixInSeg = find(SegMap==Iseg);
            SegmentsInfo(Iseg).Ind = Iseg;
            SegmentsInfo(Iseg).SegArea = length(IndPixInSeg);
            SegmentsInfo(Iseg).MeanX  = mean(MatX(IndPixInSeg));
            SegmentsInfo(Iseg).MeanY  = mean(MatY(IndPixInSeg));
            SegmentsInfo(Iseg).MinX   = min(MatX(IndPixInSeg));
            SegmentsInfo(Iseg).MaxX   = max(MatX(IndPixInSeg));
            SegmentsInfo(Iseg).MinY   = min(MatY(IndPixInSeg));
            SegmentsInfo(Iseg).MaxY   = max(MatY(IndPixInSeg));
            SegmentsInfo(Iseg).RangeX = max(MatX(IndPixInSeg)) - min(MatX(IndPixInSeg));
            SegmentsInfo(Iseg).RangeY = max(MatY(IndPixInSeg)) - min(MatY(IndPixInSeg));

            % fit line to segment
            %Par = polyfit(MatX(IndPixInSeg),MatY(IndPixInSeg),1)
            YY = MatY(IndPixInSeg);
            Err= FlatMat(IndPixInSeg);   % Scaling the errors according to the Flat value
            Err(Err<0) = 0;              % This is required as some flat pixels near the edges maybe negative
            XX = MatX(IndPixInSeg) - SegmentsInfo(Iseg).MeanX;
            NN = length(IndPixInSeg);
            H1  = [       XX, ones(NN,1)];
            H2  = [XX.^2, XX, ones(NN,1)];

            [Par1,ParErr1] = lscov(H1,YY,Err);
            [Par2,ParErr2] = lscov(H2,YY,Err);
            SegmentsInfo(Iseg).Par1    = Par1;
            SegmentsInfo(Iseg).ParErr1 = ParErr1;        
            SegmentsInfo(Iseg).Par2    = Par2;
            SegmentsInfo(Iseg).ParErr2 = ParErr2;


            % verify that Par1 and Par2 solutions are not too different
            FlagGood = abs(Par2(1))<abs(Par1(1).*0.1) && abs(Par1(1)-Par2(2))<0.01 && abs(Par1(2)-Par2(3))<0.2;
            SegmentsInfo(Iseg).FlagGood = FlagGood;

            % plot (debuging)
            %ds9_disp(InPar.IFU_FlatImageName);
            Step = 30;
            Xplot = (SegmentsInfo(Iseg).MinX-InPar.XstartOffset:Step:SegmentsInfo(Iseg).MaxX+InPar.XendOffset).'-SegmentsInfo(Iseg).MeanX;
            Yplot1 = polyval(SegmentsInfo(Iseg).Par1,Xplot);
            %ds9_plotregion(Xplot+SegmentsInfo(Iseg).MeanX,Yplot1,'Type','circle','Size',1,'color','green');
            Yplot2 = polyval(SegmentsInfo(Iseg).Par2,Xplot);
            %ds9_plotregion(Xplot+SegmentsInfo(Iseg).MeanX,Yplot2,'Type','circle','Size',1,'color','red');

            FlagGood = isempty(find(abs(Yplot1-Yplot2)>2,1));
            SegmentsInfo(Iseg).FlagGood = FlagGood;
        end
        %if (InPar.MultiThred>1),
        %   matlabpool('close');
        %end
        %toc

        % fit a surface to Slope parameter
        Xa = [SegmentsInfo.MeanX].';
        Ya = [SegmentsInfo.MeanY].';
        AllPar1 = [SegmentsInfo.Par1].';
        MeasuredSlope = AllPar1(:,1);

        [Res,FunStr]=fit_2d_polysurface(Xa,Ya,MeasuredSlope,[],'X',InPar.PolyXslope,...
                                                               'Y',InPar.PolyYslope,...
                                                               'XY',InPar.PolyXYslope);
        SlopeSurface.Res    = Res;
        SlopeSurface.FunStr = FunStr;
        
        %ParSlopeSurface = Res.Par;
        PredictedSlope  = Res.PredZ;
        %Resid            = Res.Resid;
        PredSlopeCell    = num2cell(PredictedSlope);
        [SegmentsInfo(1:1:Nseg).PredictedSlope] = deal(PredSlopeCell{:});
        
        
        % identify possible bad segments
        [SegmentsInfo(1:1:Nseg).BadSegmentFlag] = deal(num2cell([SegmentsInfo.RangeX]>170));
    otherwise
        error('Unknown FitSlope option');
end

%switch lower(InPar.SlopeType)
%    case 'measured'
%        % AllPar1 is already defined
%    case 'smooth'
%        AllPar1(:,1) = [SegmentsInfo(1:1:Nseg).PredictedSlope].';
%    otherwise
%        error('Unknown SlopeType option');
%end

%find([SegmentsInfo.FlagGood]==0)
    
