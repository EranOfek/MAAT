function [BiasSim,Summary,IsBias]=bias(Sim,varargin)
% Generate a super bias image
% Package: @SIM
% Description: Given a list of SIM images, look for the bias images,
%              and calculate the super bias image.
% Input  : - SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'IsBias'    - Vector of logical flags, or vector of indices,
%                          indicating, for each image in the SIM array,
%                          if its a bias image.
%                          If empty then attempt to find the bias images
%                          using the SIM/isbias.m function.
%                          If value is (a scalar) true then assume all
%                          the input images are bias images.
%                          Default is empty.
%            'IsBiasPar' - Cell array of additional arguments to pass to
%                          isbias.m. See isbias.m for options.
%                          Default is {}.
%            'ExecField' - Cell array of SIM objects fields which to caodd.
%                          Default is {'Im'}.
%            'CoaddMethod'- Method to coadd the bias images to a superbias.
%                          See SIM/coadd.m for options.
%                          Default is 'meanclip'.
%            'CoaddFunPar'- Cell array of additional arguments to pass
%                          to the coaddition function.
%                          Default is
%                          {'MeanFun',@nanmean,'RejectMethod','std','Reject',[4 4],'MaxIter',1}.
%            'StdMethod'  - Std calculation method. Default is @nanstd.
%            'EstimatorMean' - Function handle to use in the mean bias
%                          level estimation. Default is @nanmedian.
%            'EstimatorRN' - Function handle to use in the mean eadnoise
%                          level estimation. Default is @nanmedian.
%            'Noise0PixPar' - Noise0 pixels definition. Dead pixels (i.e.,
%                          Noise0 pixels) are defined to have std which is
%                          lower than the mean std bias level multiplied
%                          by this parameter. Default is 0.03.
%            'NoisyPixPar' - Noisy pixels definition. Noisy pixels are
%                          defined to have std which is higher than the
%                          mean std bias level multiplied by this
%                          parameter. Default is 3.
%            'AnomPixPar' - Anomalous pixels definition [low, high].
%                          Anomalous value pixels are defined to be
%                          pixels which value is lower or higher than the
%                          mean bias level +/- the parameter multiplied by
%                          readnoise in ADU units.
%                          Default is [3.5 3.5].
%            'FlarePixPar' - Flaring pixel definition. Number of sigmas
%                          above std.
%                          A flaring pixels is defined to have max/std
%                          larger than this parameter.
%                          Default is 10.
%            'BitType'   - Mask bit type (if not defined).
%                          Default is 'uint32'.
%            'MaskDic'   - Mask dictionary function handle to use.
%                          Default is @MASK.def_bitmask_pipeline.
%            'Bit_CombFun'- Function to use in the bit combining operation.
%                          Default is @bitor.
%            'Bit_BiasNoise0' - The mask bit name in which to indicate
%                          Noise0 (dead; low std) pixels.
%                          Default is 'Bit_BiasNoise0'.
%            'Bit_BiasNoisy' - The mask bit name in which to indicate
%                          Noisy (large std) pixels.
%                          Default is 'Bit_BiasNoisy'.
%            'Bit_BiasAnom' - The mask bit name in which to indicate
%                          anomalous value (low/high value) pixels.
%                          Default is 'Bit_BiasAnom'.
%            'Bit_BiasFlare' - The mask bit name in which to indicate
%                          flaring pixels.
%                          Default is 'Bit_BiasFlare'.
%            'PopHeader' - Populated basic header. Default is true.
%            'ReplaceKey'- A three column cell array of {key,val,comment},
%                          or an HEAD object which keywords to replace
%                          or add to the SIM header.
%                          Default is {}.
%            'AddKey'    - Like 'ReplaceKey', but adding keywords,
%                          without replacment.
%                          Default is {}.
%            'CoaddHeader'-Populate header with keywords related to
%                          the coaddition process. See SIM/coadd.m for
%                          details. Default is true.
% Output : - A SIM object containing the super bias image. The bias std
%            image is stored in the 'ErrIm' field, while the 'WeightIm'
%            field is populated with an image counting the number of good
%            values used for the mean estimation in each pixel.
%            The 'Mask' field contains the bit mask image.
%          - Structure containing additional meta data. Avaialble fields
%            are:
%            'MeanBias' - Mean bias level.
%            'RN_ADU'   - Mean std level (i.e., Readnoise in ADU units).
%            'BiasNoise0Map'- Image of logoicals indicating Noise0 pixels
%                             (true) or good pixels (false).
%            'BiasNoisyMap' - Image of logoicals indicating Noisy pixels
%                             (true) or good pixels (false).
%            'BiasAnomMap'  - Image of logoicals indicating Anomalous
%                             pixels (true) or good pixels (false).
%          - Vector of logicals (or indices) indicating bias images among
%            the input images.
% See also: isbias.m, debias.m, bias_overscan.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [BiasSim,Summary]=bias(Sim);
%          [BiasSim,Sum]=bias(S,'IsBias',F,'CoaddMethod',@nanmean','CoaddFunPar',{});
% Reliable: 2
%--------------------------------------------------------------------------

ImageField     = 'Im';
ErrField       = 'ErrIm';
MaskDicField   = 'MaskDic';
ImageNameField = 'ImageFileName';

DefV.IsBias             = [];
DefV.IsBiasPar          = {};
DefV.ExecField          = {ImageField};
DefV.CoaddMethod        = 'meanclip';
DefV.CoaddFunPar        = {'MeanFun',@nanmean,'RejectMethod','std','Reject',[4 4],'MaxIter',1};
DefV.EstimatorMean      = @nanmedian;    % function to estimate global meab bias level
DefV.EstimatorRN        = @nanmedian;    % function to estimate RN from bias/std image
DefV.StdMethod          = @nanstd;
DefV.Noise0PixPar       = 0.03;
DefV.NoisyPixPar        = 3;
DefV.AnomPixPar         = [3.5 3.5];
DefV.FlarePixPar        = 10;
DefV.BitType            = 'uint32';
DefV.MaskDic            = @MASK.def_bitmask_pipeline;
DefV.Bit_CombFun        = @bitor;
DefV.Bit_BiasNoise0     = 'Bit_BiasNoise0';
DefV.Bit_BiasNoisy      = 'Bit_BiasNoisy';
DefV.Bit_BiasAnom       = 'Bit_BiasAnom';
DefV.Bit_BiasFlare      = 'Bit_BiasFlare';
DefV.OutBiasName        = 'Bias.fits';
DefV.PopHeader          = true;
DefV.ReplaceKey         = {};
DefV.AddKey             = {};
DefV.CoaddHeader        = true;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% set bias images to use for superbias construction
if (isempty(InPar.IsBias))
    % Search for bias images
    IsBias = isbias(Sim,InPar.IsBiasPar{:});
elseif (islogical(InPar.IsBias))
    if (numel(InPar.IsBias)==1)
        % Assume all input images are bias images
        IsBias = true(size(Sim));
    else
        % Assume user provided flags indicating which image is a bias image
        IsBias = InPar.IsBias;
    end
else
    % Assume user provided a list of indices of bias images
    IsBias = InPar.IsBias;
end
    
if (all(~IsBias))
    % No bias images found
    error('No bias images found');
end

% construct the super bias image
% zero and scale are set to none
[BiasSim] = coadd(Sim(IsBias),'ExecField',InPar.ExecField,...
                              'CoaddMethod',InPar.CoaddMethod,...
                              'CoaddFunPar',InPar.CoaddFunPar,...
                              'MethodZero','none',...
                              'MethodScale','none',...
                              'StdMethod',InPar.StdMethod,...
                              'PopHeader',InPar.PopHeader,...
                              'ReplaceKey',InPar.ReplaceKey,...
                              'AddKey',InPar.AddKey,...
                              'CoaddHeader',InPar.CoaddHeader);


                          
% Indicate to which filter this flat corresponds
% Note that if OutFlatName contains %s than it will be replaced by the Filter name
BiasSim.(ImageNameField) = sprintf(InPar.OutBiasName);                                     

% summary and statistics
% estimate mean bias level
Summary.MeanBias  = InPar.EstimatorMean(BiasSim.(ImageField)(:));
% estimate readnoise in ADU units
Summary.RN_ADU    = InPar.EstimatorRN(BiasSim.(ErrField)(:));



%---------------------------
%--- look for bad pixels ---
%---------------------------

% 1. look for dead pixels - defined by very low std.
if (~isempty(InPar.Bit_BiasNoise0) || nargout>1)
    Summary.BiasNoise0Map  = BiasSim.(ErrField) < (Summary.RN_ADU.*InPar.Noise0PixPar);
end

% 2. look for noisy pixels - defined by very high std
if (~isempty(InPar.Bit_BiasNoisy) || nargout>1)
    Summary.BiasNoisyMap   = BiasSim.(ErrField) > (Summary.RN_ADU.*InPar.NoisyPixPar);
end

% 3. look for anomalous bias level - defined by value very different than
% mean value
if (~isempty(InPar.Bit_BiasAnom) || nargout>1)
    Summary.BiasAnomMap    = BiasSim.(ImageField) < (Summary.MeanBias - InPar.AnomPixPar(1).*Summary.RN_ADU ) | ...
                             BiasSim.(ImageField) > (Summary.MeanBias + InPar.AnomPixPar(2).*Summary.RN_ADU);
end

% 4. look for flaring pixels
% defined by normal rstd, but some anamalous rare values
if (~isempty(InPar.Bit_BiasFlare) || nargout>1)
    [MaxSim] = coadd(Sim(IsBias),'ExecField',InPar.ExecField,...
                              'CoaddMethod','max',...
                              'MethodZero','none',...
                              'MethodScale','none');
                             
    Summary.BiasFlareValue = MaxSim.(ImageField)./BiasSim.(ErrField);
    Summary.BiasFlareMap   = Summary.BiasFlareValue>InPar.FlarePixPar;
end



if (isempty(BiasSim.(MaskDicField)))
    BiasSim.(MaskDicField) = InPar.MaskDic;
end
          
% Construct the mask image
if (~isempty(InPar.Bit_BiasNoise0))
    BiasSim = bitmask_set(BiasSim,Summary.BiasNoise0Map,InPar.Bit_BiasNoise0,InPar.BitType,InPar.Bit_CombFun);
end
if (~isempty(InPar.Bit_BiasNoisy))
    BiasSim = bitmask_set(BiasSim,Summary.BiasNoisyMap, InPar.Bit_BiasNoisy, InPar.BitType,InPar.Bit_CombFun);
end
if (~isempty(InPar.Bit_BiasAnom))
    BiasSim = bitmask_set(BiasSim,Summary.BiasAnomMap,  InPar.Bit_BiasAnom,  InPar.BitType,InPar.Bit_CombFun);
end



    
    
