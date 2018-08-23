function [Sim,Back]=spec_skysub(ExtractedSpec,varargin)
%--------------------------------------------------------------------------
% spec_skysub function                                              ImSpec
% Description: Subtract the sky from a 2-dimensional spectrum.
% Input  : - The extraected 2-D spectrum in which the trace is aligned
%            with the X-axis and the spatial direction is the Y-axis.
%            This can be a single image in one of the following forms:
%            (1) A structure array (SIM).
%            (2) A file name.
%            (3) A matrix.
%            (4) A cell array with a single file name string.
%            (5) A cell array with a single matrix image.
%            See image2sim.m for options.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ReplaceIm'    - In the output image, replace the image
%                             with the background subtracted image
%                             {false|true}. Default is true.
%            'TraceCenter' - The pixel index of the trace center in
%                            the extracted spectrum. Default is empty.
%                            If empty, then this is set to be the centeral
%                            pixel of the extraction.
%            'BackRegion'  - Regions from which to extract background.
%                            These are indices range of spatial pixels
%                            in the extracted spectrum.
%                            Usually, this will have two lines, and two
%                            columns. The first line for the left region,
%                            and the second line for the right region.
%                            First column for left-hand side of a region,
%                            and second column for the right-hand side.
%                            Default is empty.
%            'FracBackRange'-If background region is not provided, than
%                            the specified fraction of left-hand
%                            and right hand pixels will be used.
%                            Default is 0.2.
%            'BackAlgo'    - Onr of the follwoing background estimation
%                            algorithm:
%                            'median' - Median of each dispersion-axis.
%                            'mean'   - Mean of each dispersion-axis.
%                            'poly'   - Polynomial fitting for each
%                                       dispersion position.
%                            'polysc' - Polynomial fitting for each
%                                       dispersion position with clipping
%                                       (default).
%            'BackPolyDeg'  - Vector of degrees of the background
%                             polynomial. Default is [0 1 2]
%                             (i.e., second order polynomial).
%            'ClipMethod'   - Clipping method for background.
%                             Default is 'MinMax'. See clip_resid.m for
%                             options.
%            'Clip'         - Clip vector for clip_resid.m.
%                             Default is [1 3].
%            'Gain'         - CCD gain. Default is 1.
%            'RN'           - CCD readnoise. Default is 10 e-.
%            'Int'          - Interactive mode {false|true}.
%                             Default is false.
%            'Range'        - Range of dispersion axis which to
%                             collapse and dispaly in the interactive mode.
%                             Default is empty. If empty, use all.
%            'Ncollapse'    - Number of collapsed regions along the
%                             dispersion axis. Default is 3.
%            'AdjPlotBack'  - Adjust background of the collapsed regions
%                             by one of the following methods:
%                             {'div'|'sub'|'no'}. Default is 'no'.
%            'CollapseAlgo' - Collapse algorithm:
%                             'optimal' - Mean weighted by inverse variance
%                                         (including readnoise). Default.
%                             'mean'    - Mean.
%                             'median'  - Median.
%                             'quantile'- Quantile with probability given by
%                                         Prct.
%            'Prct'         - Fraction for the quantile calculation for
%                             spec_collapse_dispaxis.m.
%                             Default is 0.9. This may be good for selecting
%                             sources with strong emission lines.
% Output : - Structre array containing the following fields:
%            .Im       - Background subtracted image
%            .BackIm   - Background image
%            .ErrIm    - Error image
%          - Structure with background information.
%            .BackRMS  - The rms of background best fit per dispersion
%                        position.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                     Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [BackSubSpec,BackInfo]=spec_skysub(ExtractedSpec,'BackRegion',[20 30; 70 80]);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.array.*

ImageField  = 'Im';
%HeaderField = 'Header';
%FileField   = 'ImageFileName';
%MaskField   = 'Mask';
BackImField = 'BackIm';
ErrImField  = 'ErrIm';

%InFileField = 'InputFileName';

DefV.ReplaceIm     = true;       % replace output image with background subtracted image
DefV.TraceCenter   = [];
DefV.BackRegion    = [];
DefV.FracBackRange = 0.2;         % default fraction of left and right back range in units of the extracted spectrum width
DefV.BackAlgo      = 'polysc';    % background algorithm {'median'|'mean'|'poly'}
DefV.BackPolyDeg   = [0 1 2];       
DefV.Gain          = 1;
DefV.RN            = 10;
% interactive mode
DefV.Int           = false;
DefV.Range         = [];          % range to collapse
DefV.Ncollapse     = 3;           % number of collpases in range 
DefV.AdjPlotBack   = 'div';       % adjust plot background {'div'|'sub'|'no'}
% spec_collapse_dispaxis.m parameters
DefV.CollapseAlgo  = 'optimal';   % {'optimal'|'mean'|'median','quantile'}
DefV.Prct          = 0.9;
DefV.ClipMethod    = 'MinMax';
DefV.Clip          = [1 3];
% mask
%DefV.Bit_BackRMS   = @def_bitmask_specpipeline; %5;
%DefV.BackRMS_ExFac = 2;
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

%InPar.Bit_BackRMS = get_bitmask_def(InPar.Bit_BackRMS,'Bit_BackRMS');
%InPar.BackRMS_ExFac = get_bitmask_def(InPar.BackRMS_ExFac,'BackRMS_ExFac');

Sim = image2sim(ExtractedSpec);
% if (isnumeric(ExtractedSpec)),
%     Sim.(ImageField) = ExtractedSpec;
% elseif (isstruct(ExtractedSpec)),
%     Sim = ExtractedSpec;
% else
%     error('Unknown ExtractedSpec type option');
% end

Size = size(Sim.(ImageField));
if (isempty(InPar.Range))
    InPar.Range = [1 Size(2)];
end
    
if (isempty(InPar.TraceCenter))
    % if TraceCenter is not provided than assume its in the center
    % of the extracted region
    InPar.TraceCenter = (Size(1) + 1)./2;
end
%VecSpat  = (-(Size(1)-1).*0.5:1:(Size(1)-1).*0.5).';
VecSpatI = (1:1:Size(1)).';

if (isempty(InPar.BackRegion))
    % set the default background range using the FracBackRange parameter
    InPar.BackRegion      = zeros(2,2);
    InPar.BackRegion(1,:) = [1 floor(Size(1).*InPar.FracBackRange)];
    InPar.BackRegion(2,:) = [ceil(Size(1).*(1-InPar.FracBackRange)) Size(1)];
end

if (InPar.Int)
    %--- select background interactively ---
    % select collapse ranges
    ColSize = range(InPar.Range)./InPar.Ncollapse;
    InPar.Range = round([[InPar.Range(1):ColSize:(InPar.Range(2)-ColSize)].',[InPar.Range(1)+ColSize:ColSize:InPar.Range(2)].']);
    
    Colors = generate_colors(InPar.Ncollapse);
    String = cell(1,InPar.Ncollapse);
    for Ic=1:1:InPar.Ncollapse
       Collapse = spec_collapse_dispaxis(Sim.(ImageField),varargin{:},'DispDir','x','Range',InPar.Range(Ic,:));
       switch lower(InPar.AdjPlotBack)
           case 'div'
               PlotVector = Collapse{1}./nanmedian(Collapse{1});
           case 'sub'
               PlotVector = Collapse{1} - nanmedian(Collapse{1});
           case {'no','none'}
               PlotVector = Collapse{1};
           otherwise
               error('Unknown AdjPlotBack option');
       end
       H=plot(VecSpatI,PlotVector); %, '-','Color',Colors(Ic,:))
       set(H,'Color',Colors(Ic,:));
       hold on;
       String{Ic} = sprintf('Disp %d..%d',InPar.Range(Ic,:));
    end
    legend(String{:});
    H = xlabel('Spatial position [pix]');
    set(H,'FontSize',18);
    H = ylabel('Normalized counts');
    set(H,'FontSize',18);

    fprintf('\n\n--- spec_skysub.m interactive mode ---\n\n');
    fprintf('Select Background interactively\n');
    fprintf('Select left background region\n');
    RectL = getrect;
    fprintf('Select right background region\n');
    RectR = getrect;
    InPar.BackRegion = [RectL(1), RectL(1)+RectL(3); RectR(1), RectR(1)+RectR(3)];
    
end

% measure background

IndBack   = Util.array.find_ranges(VecSpatI,InPar.BackRegion);

switch lower(InPar.BackAlgo)
    case 'median'
        BackImage = repmat(nanmedian(Sim.(ImageField)(IndBack,:)),Size(1),1);
        BackRMS   = nanstd(Sim.(ImageField)(IndBack,:));
    case 'mean'
        BackImage = repmat(nanmean(Sim.(ImageField)(IndBack,:)),Size(1),1);
        BackRMS   = nanstd(Sim.(ImageField)(IndBack,:));
    case 'poly'
        
        % construct design matrix
        Y     = Sim.(ImageField)(IndBack,:);
        Npar  = length(InPar.BackPolyDeg);
        Neq   = length(IndBack);
        H     = zeros(Neq,Npar);
        Hfull = zeros(Size(1),Npar);
        for Ipar=1:1:Npar,
           H(:,Ipar)     = IndBack.'.^InPar.BackPolyDeg(Ipar);
           Hfull(:,Ipar) = VecSpatI.'.^InPar.BackPolyDeg(Ipar);
        end
        % fit
        Par = H\Y;
        
        BackImageB = H*Par;
        ResidB     = Sim.(ImageField)(IndBack,:) - BackImageB;
        BackRMS    = nanstd(ResidB);
        BackImage  = Hfull*Par;
        
    case 'polysc'
        
        % construct design matrix
        Npar  = length(InPar.BackPolyDeg);
        Neq   = length(IndBack);
        H     = zeros(Neq,Npar);
        Hfull = zeros(Size(1),Npar);
        for Ipar=1:1:Npar,
           H(:,Ipar)     = IndBack.'.^InPar.BackPolyDeg(Ipar);
           Hfull(:,Ipar) = VecSpatI.'.^InPar.BackPolyDeg(Ipar);
        end
        
        % run over all dispersions
        Par = zeros(Npar,Size(2));
        for Idisp=1:1:Size(2),
            Y = Sim.(ImageField)(IndBack,Idisp);
            Par(:,Idisp) = H\Y;
            Resid = Y - H*Par(:,Idisp);
            [~,~,Flag] = clip_resid(Resid,'Method',InPar.ClipMethod,'Clip',InPar.Clip);
            Par(:,Idisp) = H(Flag,:)\Y(Flag);
            
        end
        BackImageB = H*Par;
        ResidB     = Sim.(ImageField)(IndBack,:) - BackImageB;
        BackRMS    = nanstd(ResidB);
        BackImage  = Hfull*Par;
        
    otherwise
        error('Unknown BackAlgo option');
end


% background subtracted image
ImageBS    = Sim.(ImageField) - BackImage;   % background subtracted image

% theoretical noise image
ImageNoise = sqrt(Sim.(ImageField)./InPar.Gain + InPar.RN.^2);

% Background RMS in:
%BackRMS
%MedianRMS = median_sigclip(BackRMS);

% Flag of back rms larger than back rms threshold
%Flag_BackRMS = BackRMS>MedianRMS.*InPar.BackRMS_ExFac;
%Flag_BackRMS = repmat(Flag_BackRMS,Size(1),1);
%
%if (isfield(Sim,MaskField)),
%    Sim.(MaskField)=maskflag_set(Sim.(MaskField),[],InPar.Bit_BackRMS,Flag_BackRMS);
%else
%    Sim.(MaskField)=maskflag_set([],[],InPar.Bit_BackRMS,Flag_BackRMS);
%end

Sim.(BackImField) = BackImage;
Sim.(ErrImField)  = ImageNoise;

if (InPar.ReplaceIm)
    Sim.(ImageField) = ImageBS;
end
Back.RMS = BackRMS;

% 
% 
% 
% clf
% imshow(ImageNoise,[0 30]);
% jj
% 
% figure(1)
% plot(BackRMS)
% figure(2);
% plot(BackImage(51,:))
% jj
% to add:
% calculate theoretical noise using Gain and RN and sqrt(N)
% calculate RMS noise
% propogate mask?
% bitset noise sky
% bitset skylines?

%plot(BackImage(51,:))
% clf
% plot(ImageBS(35,:))
% hold on
% plot(sqrt(BackImage(51,:)),'r-')
% clf
% imshow(ImageBS,[-30 100])
% fitswrite(ImageBS,'try2.fits')
% 


    
    
    
    
