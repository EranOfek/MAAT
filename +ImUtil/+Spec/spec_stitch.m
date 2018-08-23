function [MList]=spec_stitch(List1,List2,varargin)
%--------------------------------------------------------------------------
% spec_stitch function                                              ImSpec
% Description: Given two spectra of a single object taken at different
%              wavelength (e.g., blue side and red side), with or without
%              an overlap between the spectra, this program find an
%              optimal flux ratio shift between the two spectra,
%              and return the stitched spectrum.
% Input  : - A matrix containing a blue-side spectrum to stitch
%            [Wave, Flux, Err,...]. The wavelength should be sorted.
%            The matrix may containin an arbitrary number of columns
%            (e.g., background). This columns will be interpolated
%            and corrected as well.
%          - A matrix containing a red-side spectrum to stitch
%            [Wave, Flux, Err,...]. The wavelength should be sorted.
%          * Arbitrary number of ...,key,val,... arguments.
%            The following keywords are available:
%            'Algo'   - One of the following stitching algorithms:
%                       'lsq'  - Find a shift based on linear least
%                                sqaure fitting of the overlapping
%                                sections (default).
%                       'poly' - Find a shift that minimize the rms
%                                of an high order polynomial fitted
%                                to both spectra.
%            'PolyOrder' - Order of polynomial in the 'poly' algorithm.
%                       Default is 4.
%            'FluxRatioVec' - Vector of flux ratios to test in the
%                       'poly' algorithm. Default is
%                       1./logspace(-1,1,30).'.
%            'Merge'  - Merging method for data. Options are:
%                       'mean' - use the mean of the blue and red overlap.
%                       'vmean'- variance weighted mean. Default.
%                       'smean'- std weighted mean.
%                       '1'  - use blue side spectrum.
%                       '2'  - use red side spectrum.
%            'IntList'- Which input list to interpolate {1|2}. Default
%                       is 1.
%            'BM1'    - A vector of bit masks for the first spectrum.
%                       Default is empty.
%            'BM2'    - A vector of bit masks for the second spectrum.
%                       Default is empty.
%            'BMb'    - Bit mask merging behavior {'and'|'or'}.
%                       Default is 'or'.
%            'AddBM'  - The index of a bit to add in the overlapping
%                       region. If empty matrix then do not add a bit.
%                       Default is empty.
%            'TypeBM' - Type of bit mask {'uint8'|'uint16'|'uint32'}.
%                       Default is 'uint16'.
%            'InterpMethod' - Interpolation method for spectra.
%                       See interp1.m for options. Default is 'linear'.
%            'ClipPar'- Cell array of parameters to pass to clip_resid.m
%                       for sigma clipping. See clip_resid for options
%                       and defaults. Default is no sigma clipping.
%            'MaxBlue'- Maximum wavelength to use in the blue side.
%                       Default is Inf. This can be use to truncate
%                       the blue spectrum before stitching.
%            'MinRed' - Minimum wavelength to use in the red side.
%                       Default is -Inf. This can be use to truncate
%                       the red spectrum before stitching.
% Output : - Two column matrix with the stitched spectedum.
%          - A vector of bit mask.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [ListM,BM]=spec_stitch(B,R,'AddBM',1);
% Reliable: 2
%--------------------------------------------------------------------------

Col.Wave = 1;
Col.Flux = 2;
Col.Err  = 3;


DefV.MedFiltSize    = 10;
DefV.Prctile        = 0.99;

DefV.Algo           = 'lsq'; %'lsq';     % {'lsq'|'poly'}
DefV.BM1            = [];
DefV.BM2            = [];
DefV.BMb            = 'or';
DefV.AddBM          = [];
DefV.TypeBM         = 'uint16';
DefV.Merge          = 'mean';
DefV.IntList        = 1;
DefV.InterpMethod   = 'linear';
DefV.ClipPar        = {};
DefV.MaxBlue        = Inf;
DefV.MinRed         = -Inf;
DefV.ClipPar        = {};

DefV.FluxRatioVec   = 1./logspace(-1,1,30).';
DefV.PolyOrder      = 4;
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

% extra columns:
[N1,Ncol] = size(List1);
[N2,~]    = size(List2);

if (isempty(InPar.BM1)),
    InPar.BM1 = zeros(N1,1,InPar.TypeBM);
end
if (isempty(InPar.BM2)),
    InPar.BM2 = zeros(N2,1,InPar.TypeBM);
end

% remove bad regions from spectra
List1 = List1(List1(:,Col.Wave)<InPar.MaxBlue,:);
List2 = List2(List2(:,Col.Wave)>InPar.MinRed,:);

% calculate overlap region
MinOverlapWave = min(List2(1,Col.Wave));
MaxOverlapWave = max(List1(end,Col.Wave));


if (MaxOverlapWave>MinOverlapWave),
    % there is overlap
    Overlap = true;
else
    % no overlap
    Overlap = false;
end

% interpolate list to uniform grids
StepWave = min(min(diff(List1(:,Col.Wave))),min(diff(List2(:,Col.Wave))));
MinWave  = min(min(List1(:,Col.Wave)),min(List2(:,Col.Wave)));
MaxWave  = max(max(List1(:,Col.Wave)),max(List2(:,Col.Wave)));
Wave     = (MinWave:StepWave:MaxWave).';

% uniform spectra
UList1    = [Wave, interp1(List1(:,Col.Wave),List1(:,Col.Flux),Wave,InPar.InterpMethod)];
UList2    = [Wave, interp1(List2(:,Col.Wave),List2(:,Col.Flux),Wave,InPar.InterpMethod)];

%Err1      = 


FlagOverlap = ~isnan(UList1(:,Col.Flux)) & ~isnan(UList2(:,Col.Flux));




% modify algo if no overlap
if (~Overlap)
    switch lower(InPar.Algo)
        case {'lsq'}
            warning('No overlap between two spectra switch Algo to poly');
            InPar.Algo = 'poly';
        otherwise
            % do nothing
    end
end

switch lower(InPar.Algo)
    case 'lsq'

        % generate Noise estimaate
        Noise1 = stdfilt1(UList1(:,Col.Flux),InPar.MedFiltSize);
        Noise2 = stdfilt1(UList2(:,Col.Flux),InPar.MedFiltSize);

        % remove data with high noise
        FlagGoodNoise = Noise1(FlagOverlap)<prctile(Noise1(FlagOverlap),InPar.Prctile.*100) & ...
                        Noise2(FlagOverlap)<prctile(Noise2(FlagOverlap),InPar.Prctile.*100);
        
        
        IndOverlap    = find(FlagOverlap);
        IndUse        = IndOverlap(FlagGoodNoise);
       
        Y      = UList1(IndUse,Col.Flux);
        H      = UList2(IndUse,Col.Flux);
        V      = Noise1(IndUse).^2 + Noise2(IndUse).^2;
        
        [Factor,FactorErr] = lscov(H,Y,1./V);
        
        
    otherwise
        error('Unknown Algo option');
end




switch lower(InPar.Merge)
    case 'mean'
        MList = [Wave, nanmean([UList1(:,Col.Flux),Factor.*UList2(:,Col.Flux)],2)];
      
  
    otherwise
        error('Unknown Merge option');
end

