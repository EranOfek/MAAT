function [CoaddSim,StdSim,WeightSim]=coadd(Sim,varargin)
% Simple coaddition of images in a SIM object.
% Package: @SIM
% Description: Simple (pixel by pixel) coaddition of SIM object images.
% Input  : - A multi element SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - The SIM fields on which to execute the
%                          coaddition. Default is {'Im'}.
%            'CoaddMethod'- Either a function handle or a string.
%                          The function is of the form Fun(Cube,AddPar{:})
%                          (e.g., sum(Cube,1)).
%                          Alternatively this can be one of the following
%                          strings (methods):
%                          'meanclip' (use image_mean_clip.m),...
%                          'mean', 'nanmean', 'rmean','sum', 'nansum',...
%                          'std', 'nanstd', 'median', 'nanmedian',...
%                          'quantile', 'iqr', 'var', 'nanvar',...
%                          'min', 'nanmin', 'max', 'nanmax', 'range',...
%                          'sqrtsumsq' - sqrt(sum(X^2)).
%                          'count' - count not NaN values.
%            'CoaddFunPar'- Cell array of additional arguments to pass to
%                          coadd function (if CoaddMethod is a function
%                          handle).
%            'MethodZero'- Images additive scaling method prior
%                          to coaddition. See scale.m for options.
%                          Default is 'none'.
%            'ZeroPar'   - Additional arguments to pass to scale.m
%                          for the zero subtraction. Default is {}.
%            'MethodScale'- Images multiplicative scaling method prior
%                          to coaddition. See scale.m for options.
%                          Default is 'none'.
%                          Note that zero subtraction comes prior to
%                          scaling.
%            'ScalePar'  - Additional arguments to pass to scale.m
%                          for the scale multiplication. Default is {}.
%            'ImDim'     - The image index dimension in the cube.
%                          Default is 1.
%            'StdMethod' - Method for calculating the image error (std).
%                          If empty then do not calculate std.
%                          Default is @nanstd.
%            'WeightMethod' - Method for calculating the image weight.
%                          If empty then do not calculate weight.
%                          Default is 'count'.
%            'MaskCombFun' - Method by which to coadd the mask images.
%                          Default is @bitor. See mask_combine.m for
%                          options. If empty, the  don't combine the masks.
%            'PopHeader' - Populated basic header. Default is true.
%            'ReplaceKey'- A three column cell array of {key,val,comment},
%                          or an HEAD object which keywords to replace
%                          or add to the SIM header.
%                          Default is {}.
%            'AddKey'    - Like 'ReplaceKey', but adding keywords,
%                          without replacment.
%                          Default is {}.
%            'CoaddHeader'-Populate header with keywords related to
%                          the coaddition process.
%                          Keywords include:
%                          'COADD',true,...
%                          'EXPTIME',sum(ExpTime),...
%                          'NIMAGE',numel(Images),...
%                          'JD',mean(JD),...
%                          'MIN_JD',min(JD),...
%                          'MAX_JD',max(JD),...
%                          Default is true.
% Output : - A SIM object containing the coadded images.
%            If StdMethod is not empty and nargout==1 then
%            the std image is populated in the 'ErrIm' field.
%            If WeightMethod is not empty and nargout<2 then
%            the weight image is populated in the 'WeightIm' field.
%          - A SIM object with the StD image in the 'Im' field.
%          - A SIM object with the Weight image in the 'Im' field.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Coadd with all defaults
%          C=coadd(S);
%          % sigma cliping coadd
%          C=coadd(S,'CoaddMethod','meanclip');
%          % sigma cliping coadd with min/max rejection (1 low, 5 high).
%          C=coadd(S,'CoaddMethod','meanclip','CoaddFunPar',{'RejectMethod','minmax','Reject',[1 5]});
%          % coadd with std and weight in different SIM objects.
%          [Coadd,Std,Weight]=coadd(S);
%          % only coadd image without std and weight
%          [Coadd]=coadd(S,'WeightMethod',[],'StdMethod',[]);
%          % coadd with zero subtraction and scaling
%          [Coadd]=coadd(S,'MethodZero','median','MethodScale','median');
% Reliable: 2
%--------------------------------------------------------------------------

ImageField   = 'Im';
ErrField     = 'ErrIm';
WeightField  = 'WeightIm';
%MaskField    = 'Mask';
%MaskDicField = 'MaskDic';

DefV.ExecField          = {'Im'};
DefV.CoaddMethod        = @nanmean;
DefV.CoaddFunPar        = {};
DefV.MethodZero         = 'none';
DefV.ZeroPar            = {};
DefV.MethodScale        = 'none';
DefV.ScalePar           = {};
DefV.ImDim              = 1;
DefV.StdMethod          = @nanstd;  
DefV.WeightMethod       = 'count';
DefV.MaskCombFun        = @bitor;
DefV.PopHeader          = true;
DefV.ReplaceKey         = {};
DefV.AddKey             = {};
DefV.CoaddHeader        = true;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nf = numel(InPar.ExecField);
CoaddSim = SIM;
for If=1:1:Nf
    
    % subtract some constant
    % use scale.m with ScaleFun=@minus
    Sim = scale(Sim,InPar.MethodZero,InPar.ZeroPar{:},'ExecField',InPar.ExecField{If},'ScaleFun',@minus);
    
    % Multiply by Scale
    Sim = scale(Sim,InPar.MethodScale,InPar.ScalePar{:},'ExecField',InPar.ExecField{If});
    
    % Apply filter to images
    % FFU <<<<<<<<<<<<<<<<<<
    
    % Convert SIM to Cube
    Cube = sim2cube(Sim,'ExecField',InPar.ExecField{If},'ImDim',InPar.ImDim);

    [Image,StD,Weight] = coadd_cube(Cube,InPar);    
    
    CoaddSim.(InPar.ExecField{If}) = squeeze(Image);
    
    % populate the image Std
    if (~isempty(InPar.StdMethod))
        if (isempty(StD))
            % StD was not calculated by coadd_cube
            % so calculate it...
            InParStd = InPar;
            InParStd.CoaddMethod = InPar.StdMethod;
            StD      = coadd_cube(Cube,InParStd);
        end
        
        StD = squeeze(StD);
        
        if (nargout>1)
            StdSim = SIM;
            StdSim.(ImageField) = StD;
        else
            % nargout==1
            % populate StD in ErrField:
            CoaddSim.(ErrField) = StD;
        end
    end
    
    
    % populate the image NpixUse (weight)
    if (~isempty(InPar.WeightMethod))
        if (isempty(Weight))
            % NpixUse (Weight) was not calculated by coadd_cube
            % so calculate it...
            InParW = InPar;
            InParW.CoaddMethod = InPar.WeightMethod;
            Weight   = coadd_cube(Cube,InParW);
        end
        
        Weight = squeeze(Weight);
        
        if (nargout>1)
            WeightSim = SIM;
            WeightSim.(ImageField) = Weight;
        else
            % nargout==1
            % populate Weight in WeightField:
            CoaddSim.(WeightField) = Weight;
        end
    end
    
end

% combine the imask images
if (~isempty(InPar.MaskCombFun))
    % combine mask and copy it into CoaddSim:
    CoaddSim = mask2sim(mask_combine(Sim,InPar.MaskCombFun),CoaddSim);
end

% populate basic header
if (InPar.PopHeader)
    CoaddSim = pop_basic_header(CoaddSim,'ReplaceKey',InPar.ReplaceKey','AddKey',InPar.AddKey);
end

% populated keyword related to coaddition
if (InPar.CoaddHeader)
    [JD,ExpTime] = julday(Sim,'UpdateHead',false);
    KeyH = {'COADD',true,'IMAGE BASED ON COADDITION'; ...
            'EXPTIME',nansum(ExpTime),'TOTAL EXP. TIME OF COADD IMAGE'; ...
            'NIMAGE',numel(Sim),'TOTAL NUMBER OF COADD IMAGES'; ...
            'JD',mean(JD),'MEAN JD OF COADD IMAGES'; ...
            'MIN_JD',min(JD),'MINIMUM JD OF COADD IMAGES'; ...
            'MAX_JD',min(JD),'MAXIMUM JD OF COADD IMAGES'};
    CoaddSim = replace_key(CoaddSim,KeyH);
    
end


% end of function
end

    
    function [Image,StD,NpixUse] = coadd_cube(Cube,InPar)
    %---------------------------------------------
    % Cube coaddition along one of the dimensions
    %---------------------------------------------
    StD = [];
    NpixUse = [];
    if (isa(InPar.CoaddMethod,'function_handle'))
        Image = InPar.CoaddMethod(Cube,InPar.CoaddFunPar{:});
    else    
        switch lower(InPar.CoaddMethod)
            case 'meanclip'
                % mean with sigma clipping
                [Image,StD,NpixUse] = ImUtil.Im.image_mean_clip(Cube,InPar.CoaddFunPar{:},'ImDim',InPar.ImDim);
            case 'count'
                Image = sum(~isnan(Cube),InPar.ImDim);
            case 'mean'
                Image = mean(Cube,InPar.ImDim);
            case 'nanmean'
                Image = nanmean(Cube,InPar.ImDim);
            case 'rmean'
                Image = rmean(Cube,InPar.ImDim,InPar.CoaddFunPar{:});
            case 'sum'
                Image = sum(Cube,InPar.ImDim);
            case 'nansum'
                Image = nansum(Cube,InPar.ImDim);
            case 'std'
                Image = std(Cube,0,InPar.ImDim);
            case 'nanstd'
                Image = nanstd(Cube,0,InPar.ImDim);
            case 'median'
                Image = median(Cube,InPar.ImDim);
            case 'nanmedian'
                Image = nanmedian(Cube,InPar.ImDim);
            case 'quantile'
                % CoaddFunPar should contain the probability in fraction
                Image = quantile(Cube,InPar.CoaddFunPar{:},InPar.ImDim);
            case 'iqr'
                Image = iqr(Cube,InPar.ImDim);
            case 'var'
                % CoaddFunPar should contain a weight vector (inverse
                % variance)
                Image = var(Cube,InPar.CoaddFunPar{:},InPar.ImDim);
             case 'nanvar'
                % CoaddFunPar should contain a weight vector (inverse
                % variance)
                Image = nanvar(Cube,InPar.CoaddFunPar{:},InPar.ImDim);
            case 'min'
                Image = min(Cube,[],InPar.ImDim);
            case 'nanmin'
                Image = nanmin(Cube,[],InPar.ImDim);
            case 'max'
                Image = max(Cube,[],InPar.ImDim);
            case 'nanmax'
                Image = nanmax(Cube,[],InPar.ImDim);
            case 'range'
                Image = range(Cube,InPar.ImDim);
            case 'sqrtsumsq'
                Image = sqrt(nansum(Cube.^2,InPar.ImDim));
            otherwise
                error('Unknown CoaddMethod option');
        end
    end
end

