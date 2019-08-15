function [FlatSim,Summary,IsFlat]=flat(Sim,varargin)
% Create a flat field image from images in a SIM object.
% Package: @SIM
% Description: Given a list of SIM images, construct a flat image for
%              each filter. The flat image can be constrcted from dome
%              flat, twillight flat or science images (super flat).
% Input  : - SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'IsFlat'    - Vector of logical flags, or vector of indices,
%                          indicating, for each image in the SIM array,
%                          if to use it as a flat image.
%                          If empty then attempt to find the flat images
%                          using the SIM/isflat.m function.
%                          If value is (a scalar) true then assume all
%                          the input images are bias images.
%                          Default is empty.
%            'IsFlatPar' - Cell array of additional arguments to pass to
%                          isflat.m. See isbias.m for options.
%                          Default is {}.
%            'FlagMap'   - A User supplied Flag Map, for each input image,
%                          which is either a
%                          SIM in which the Im field contains logical
%                          flags, a matrix of logical flags, or a matrix
%                          of bit masks. If empty, then will attempt to use
%                          the .Mask field in the input SIM.
%                          Default is empty.
%                          In order to skip this option set 'BitMask' to 0.
%            'BitMask'   - If any of the bit mask are satisfied by the
%                          optional bit mask provided in 'FlagMap'
%                          (or the .Mask field) then 
%                          the pixel will be set to NaN so it will not be
%                          used for the flat field.
%                          E.g., this can be used to ignore stars and bad
%                          pixels to participate in the flat calculation.
%                          Default is 0.
%            'BitType'   - If the 'BitMask' in the SIM is numeric
%                          then this parameter is used to specify if the
%                          bit mask is a vector of indices ('index') or a
%                          bit mask value ('value'). Default is 'value'.
%            'ExecField' - Cell array of SIM objects fields which to coadd.
%                          Default is {'Im'}.
%            'FilterKey' - A filter keyword string or a cell array of
%                          filter keyword strings by which to group the
%                          input images. A flat field image will be
%                          constructed for each group.
%                          Default is 'FILTER'.
%            'CoaddMethod'- Method to coadd the flat images to a superflat.
%                          See SIM/coadd.m for options.
%                          Default is 'meanclip'.
%            'CoaddFunPar'- Cell array of additional arguments to pass
%                          to the coaddition function.
%                          Default is
%                          {'MeanFun',@nanmean,'RejectMethod','std','Reject',[3 3],'MaxIter',1}.
%            'StdMethod' - Method for calculating the flat image std.
%                          If empty then do not calculate std.
%                          See SIM/coadd.m for details.
%                          The 'ErrImType' argument specify if to divide
%                          this by sqrt(Nimages).
%                          Default is @nanstd.
%            'ErrImType'  - Specify the content of the 'ErrIm' field in
%                          Flat image. One of the following options:
%                          'std' - The error image calculation is defined
%                                  by the 'StdMethod'.
%                          'err' - The error image resulted from the
%                                  'StdMethod' option is divided by the
%                                  number of images used for each pixel
%                                  (i.e. weight image). 
%                          'relerr' - The error image resulted from the
%                                  'StdMethod' option is divided by the
%                                  number of images used for each pixel
%                                  (i.e. weight image) and than divided by
%                                  the value of the flat in each position
%                                  (i.e., this is the relative error).
%                          Default is 'relerr'.
%            'MethodZero' - Subtract zero level method prior to coaddition.
%                          See SIM/coadd.m for options.
%                          Default is 'none'.
%            'ZeroPar'    - Additional parameters to pass to the zero
%                          method prior to coaddition.
%                          See SIM/coadd.m for options.
%                          Default is {}.
%            'MethodScale'- Scale method prior to coaddition.
%                          See SIM/coadd.m for options.
%                          Default is '1/mean'.
%            'ScalePar'   - Additional parameters to pass to the scale
%                          method prior to coaddition.
%                          See SIM/coadd.m for options.
%                          Default is {}.
%            'PostScale'  - Post coaddition scaling method.
%                          See SIM/scale.m for options.
%                          Default is '1/mean'.
%            'PostScalePar'- Additional parameters to pass to the scale
%                          method post to coaddition.
%                          See SIM/coadd.m for options.
%                          Default is {}.
%            'OutFlatName'- The name of the flat image in the ImageFileName
%                          field in the SIM object. If the file name
%                          contains an %s, then this will be replaced with
%                          the merged filter name, where merged filter name
%                          is the string containing all the filter names
%                          from the 'FilterKey' keywords.
%                          Default is 'Flat_%s.fits'.
%            'EstimatorMean' - Function handle to use in the mean flat
%                          level, std and weight estimation reported in the
%                          Summary output. Default is @nanmedian.
%            'BitClass' - Class of the newly formed MASK image.
%                         Default is 'uint32'.
%            'MaskDic'  - Default dictionary for the flat MASK image.
%                         Default is @MASK.def_bitmask_pipeline.,
%            'Bit_CombFun' - Method by which to combine the bit mask.
%                         Default is @bitor.
%            'Bit_FlatNaN' - Name of bit mask for flat pixels which are
%                         equal NaN. Default is 'Bit_FlatNaN'.
%            'Bit_FlatLowNim' - Name of bit mask for flat pixels which are
%                         the result of a low number of input pixels.
%                         The number of minimum input pixels is defined by
%                         'Par_FlatLowNim'.
%                         Default is 'Bit_FlatLowNim'.
%            'Bit_FlatHighStd' - Name of bit mask for flat pixels which
%                         have high std.
%                         The number of maximum std is defined by
%                         'Par_FlatHighStd'.
%                         Default is 'Bit_FlatHighStd'.
%            'Bit_FlatLow' - Name of bit mask for flat pixels which have
%                         low value.
%                         The low value is defined by 'Par_FlatLow'.
%                         Defaukt is 'Bit_FlatLow'.
%            'Par_FlatLowNim' - Parameter for the 'Bit_FlatLowNim' bit.
%                         Default is 4.
%            'Par_FlatHighStd' - Parameter for the 'Bit_FlatHighStd' bit.
%                         Default is 0.03.
%            'Par_FlatLow' - Parameter for the 'Bit_FlatLow' bit.
%                         Default is 0.4.
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
% Output : - A SIM object containing the super flat image. The flat std
%            or error or relative error image (see 'ErrImType' option)
%            is stored in the 'ErrIm' field, while the 'WeightIm'
%            field is populated with an image counting the number of good
%            values used for the mean estimation in each pixel.
%            The 'Mask' field contains the bit mask image.
%          - Structure containing additional meta data. Avaialble fields
%            are:
%            'FilterCell' - A cell array of all filters composing the
%                           filter group (e.g., a group may be defined by
%                           multiple filter wheels).
%            'Filter'     - The combined filter name of the group.
%            'ImInd'      - Indices of images in the input SIM that belongs
%                           to the filter group.
%            'ScaleFactor'- The scale factor applied to the final coadded
%                           image.
%            'MeanFF'     - Mean flat level (calculated using
%                          'EstimatorMean').
%            'MinFF'      - Minimum value in flat.
%            'MaxFF'      - Maximum value in flat.
%            'MeanStD'    - Mean of the 'ErrIm' field of the flat image
%                           (calculated using 'EstimatorMean').
%            'MeanNim'    - Mean number of images participated in the
%                           flat constructiom. I.e., The 'EstimatorMean'
%                           function on the .WeightIm field.
%            'FlatNaNMap' - Image of logicals indicating NaN pixels
%                           (true) or good pixels (false).
%            'FlatLowNimMap'- Image of logicals indicating low number of
%                           images used for flat construction
%                           (true) or good pixels (false).
%            'FlatHighStdMap'- Image of logicals indicating high std pixels
%                           (true) or good pixels (false).
%            'FlatLowMap' - Image of logicals indicating low value pixels
%                           (true) or good pixels (false).
%          - Vector of logicals (or indices) indicating flat images among
%            the input images.
% See also: isflat.m, deflat.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [FlatSim,Summary]=bias(Sim);
%          [FlatSim,Sum]=bias(S,'IsFlat',F,'CoaddMethod',@nanmean','CoaddFunPar',{});
% Reliable: 2
%--------------------------------------------------------------------------

ImageField     = SIM.ImageField;
ImageNameField = SIM.FileNameField;
ErrField       = SIM.ErrField;
WeightField    = SIM.WeightField;
%MaskField      = 'Mask';
MaskDicField   = 'MaskDic';


% Default input arguments                                    
DefV.IsFlat             = [];
DefV.IsFlatPar          = {};
DefV.FlagMap            = [];
DefV.BitMask            = 0;
DefV.BitType            = 'value';  % 'value'|'index' - see bitmask_find.m
DefV.ExecField          = {SIM.ImageField};
DefV.FilterKey          = 'FILTER';
DefV.MethodScale        = '1/mean';
DefV.ScalePar           = {};
DefV.MethodZero         = 'none';
DefV.ZeroPar            = {};
DefV.CoaddMethod        = 'meanclip';
DefV.CoaddFunPar        = {'RejectMethod','std','Reject',[3 3],'MaxIter',1};
DefV.StdMethod          = @nanstd;
DefV.ErrImType          = 'relerr';          % 'std'|'err'
DefV.PostScale          = '1/mean';
DefV.PostScalePar       = {};
DefV.OutFlatName        = 'Flat_%s.fits';
DefV.EstimatorMean      = @nanmedian;    % estimator for mean properties of flat image in Summary output
DefV.BitClass           = 'uint32';
DefV.MaskDic            = @MASK.def_bitmask_pipeline;
DefV.Bit_CombFun        = @bitor;
DefV.Bit_FlatNaN        = 'Bit_FlatNaN';
DefV.Bit_FlatLowNim     = 'Bit_FlatLowNim';
DefV.Bit_FlatHighStd    = 'Bit_FlatHighStd';
DefV.Bit_FlatLow        = 'Bit_FlatLow';
DefV.Par_FlatLowNim     = 4;
DefV.Par_FlatHighStd    = 0.03;
DefV.Par_FlatLow        = 0.4;
DefV.ReplaceKey         = {};
DefV.AddKey             = {'IMTYPE','FLAT'};
DefV.CoaddHeader        = true;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% divide list of images to groups by filter name
FilterGroups = find_groups(Sim,InPar.FilterKey);
Nfilter      = numel(FilterGroups);
Nsim         = numel(Sim);

FlatSim = SIM(Nfilter,1);
if (nargout>1)
    Summary = Util.struct.struct_def({'FilterCell','Filter','ImInd','ImFlag',...
                          'ScaleFactor',...
                          'MeanFF','MaxFF','MinFF','MeanStD','MeanNim',...
                          'FlatNaNMap','FlatLowNimMap','FlatHighStdMap','FlatLowMap'},Nfilter,1);
end
for Ifilter=1:1:Nfilter
    % for each filter set of images
    ImInd         = FilterGroups(Ifilter).ptr;   % indices of images
    % Flag indicating if image belong to current group
    ImFlag        = false(Nsim,1);
    ImFlag(ImInd) = true;
    
    FilterCell    = FilterGroups(Ifilter).Content;
    % Filter is a cell array of filters - the next command will concat all
    % the filters into a single filter string
    Filter     = cell2mat(FilterCell);
    
    % search for flat images 
    if (isempty(InPar.IsFlat))
        % IsFlat is empty - use isflat.m
        IsFlat = isflat(Sim(ImInd),InPar.IsFlatPar{:});
    elseif (islogical(InPar.IsFlat) && numel(InPar.IsFlat)==1)
        % IsFlat is true - so use all images
%         IsFlat = ImInd;
        IsFlat = true(size(ImInd)); % Na'ama, 2018-06-06
        
    elseif (islogical(InPar.IsFlat) && numel(InPar.IsFlat)>1)
        % IsFalt is logicals
        IsFlat = InPar.IsFlat(ImInd);
    elseif (isnumeric(InPar.IsFlat))
        % IsFlat is indices
        IsFlat = InPar.IsFlat(ImInd);
    elseif (ischar(InPar.IsFlat) || iscellstr(InPar.IsFlat))
        % IsFlat is string or a cell of strings
        IsFlat = istype(Sim(ImInd),InPar.IsFlat);
    else
        error('Unknown IsFlat option');
    end
    % IsFlat contains either logicals or indices of images to use as flat
    
    % verify that flat images exist
    if ((islogical(IsFlat) && ~any(IsFlat)) || isempty(IsFlat))
        error('Flat images were not found for filter %s',Filter);
    end
        
    % Sim(ImInd(IsFlat)) contains the images from which to construct a flat
    ImIndFlat = ImInd(IsFlat);
    
    % mask unwanted pixels
    % set unwanted pixels with NaN
    Sim(ImIndFlat) = flag2nan(Sim(ImIndFlat),'ExecField',InPar.ExecField,...
                                             'FlagMap',InPar.FlagMap,...
                                             'BitMask',InPar.BitMask,...
                                             'BitType',InPar.BitType);
                                         
   
    % noralize all the images
    % by e.g., nanmean or nanmedian
    % and coadd them
    % Note that MaskCombFun is set to [] so mask is not propogated through
    % the flat coaddition step - THe reasoning is that bad pixels can be
    % ignored due to the previous step (flag2nan.m).
    % Propogation of the image mask should be done through the image...
    FlatSim(Ifilter) = coadd(Sim(ImIndFlat),'ExecField',InPar.ExecField,...
                                          'MethodScale',InPar.MethodScale,...
                                          'ScalePar',   InPar.ScalePar,...
                                          'MethodZero', InPar.MethodZero,...
                                          'ZeroPar',    InPar.ZeroPar,...
                                          'CoaddMethod',InPar.CoaddMethod,...
                                          'CoaddFunPar',InPar.CoaddFunPar,...
                                          'StdMethod',  InPar.StdMethod,...
                                          'MaskCombFun',[],...
                                          'ReplaceKey', InPar.ReplaceKey,...
                                          'AddKey',     InPar.AddKey,...
                                          'CoaddHeader',InPar.CoaddHeader);

                      
                                      
    % Define the type of the ErrIm field ('std' or 'err')
    switch lower(InPar.ErrImType)
        case 'err'
            % divide ErrIm by sqrt(Nimages)
            % I.e., error on the flat value
            FlatSim(Ifilter).(ErrField) = FlatSim(Ifilter).(ErrField) ./ sqrt(FlatSim(Ifilter).(WeightField));
        case 'relerr'
            % divide ErrIm by sqrt(Nimages) and by value of flat
            % I.e., relative error on the flat value
            FlatSim(Ifilter).(ErrField) = FlatSim(Ifilter).(ErrField) ./ ...
                            (sqrt(FlatSim(Ifilter).(WeightField)) .*  FlatSim(Ifilter).(ImageField));
        case 'std'
            % do nothing
        otherwise
            error('Unknown ErrImType option');
    end
    
    % Post scaling                                  
    % scale final flat image (to unity)
    [FlatSim(Ifilter),ScaleFactor] = scale(FlatSim(Ifilter),InPar.PostScale,...
                                              'ExecField',InPar.ExecField,...
                                              'FunAddPar',InPar.PostScalePar);
    
    % Scale final flat err image using the same scale factor as the flat
    % image (previous step).
    FlatSim(Ifilter).(ErrField) = FlatSim(Ifilter).(ErrField)./ScaleFactor;
    
    % set the additional header keywords for the flat image
    % including: image type, filter
    FlatSim(Ifilter) = replace_key(FlatSim(Ifilter),{InPar.FilterKey,Filter,'Merged Filter name';...
                                                     'Type','Flat','Image type'});
    
    % Set the flat image name
    % Indicate to which filter this flat corresponds
    % Note that if OutFlatName contains %s than it will be replaced by the Filter name
    FlatSim(Ifilter).(ImageNameField) = sprintf(InPar.OutFlatName,Filter);                                     
                                         
    % summary and statistics
    if (nargout>1)
        Summary(Ifilter).FilterCell  = FilterCell;  % Cell array of filter names of group
        Summary(Ifilter).Filter      = Filter;      % The combined filter name
        Summary(Ifilter).ImInd       = ImInd;       % indices of images taken with this filter
        Summary(Ifilter).ImFlag      = ImFlag;      % Flag indicating of image belong to current group
        Summary(Ifilter).ScaleFactor = ScaleFactor;
        Summary(Ifilter).MeanFF      = InPar.EstimatorMean(FlatSim(Ifilter).(ImageField)(:));
        Summary(Ifilter).MaxFF       = nanmax(FlatSim(Ifilter).(ImageField)(:));
        Summary(Ifilter).MinFF       = nanmin(FlatSim(Ifilter).(ImageField)(:));    
        Summary(Ifilter).MeanStD     = InPar.EstimatorMean(FlatSim(Ifilter).(ErrField)(:));
        Summary(Ifilter).MeanNim     = InPar.EstimatorMean(FlatSim(Ifilter).(WeightField)(:));
    end

    %---------------------------
    %--- look for bad pixels ---
    %---------------------------

    % 1. look for NaN
    % Bit_FlatNaN
    if (~isempty(InPar.Bit_FlatNaN) || nargout>1)
        Summary(Ifilter).FlatNaNMap = isnan(FlatSim(Ifilter).(ImageField));
    end
    % 2. look for pixels with low number of images
    % Bit_FlatLowNim
    if (~isempty(InPar.Bit_FlatLowNim) || nargout>1)
        Summary(Ifilter).FlatLowNimMap = FlatSim(Ifilter).(WeightField)<InPar.Par_FlatLowNim;
    end
    % 3. look for pixels with large std
    % Bit_FlatHighStd
    if (~isempty(InPar.Bit_FlatHighStd) || nargout>1)
        Summary(Ifilter).FlatHighStdMap = FlatSim(Ifilter).(ErrField)>InPar.Par_FlatHighStd;
    end
    % 4. look for low flat vlaues
    % Bit_FlatLow
    if (~isempty(InPar.Bit_FlatLow) || nargout>1)
        Summary(Ifilter).FlatLowMap = FlatSim(Ifilter).(ImageField)<InPar.Par_FlatLow;
    end
 
    % Populate the MaskDic field
%     if (isempty(FlatSim.(MaskDicField)))
%        FlatSim.(MaskDicField) = InPar.MaskDic;
%     end
    if (isempty(FlatSim(Ifilter).(MaskDicField))) % Na'ama, 2018-06-06
       FlatSim(Ifilter).(MaskDicField) = InPar.MaskDic;
    end
    
    % Construct the mask image
    if (~isempty(InPar.Bit_FlatNaN))
        %FlatSim = bitmask_set(FlatSim,Summary.FlatNaNMap,InPar.Bit_FlatNaN,InPar.BitClass,InPar.Bit_CombFun);
        FlatSim(Ifilter) = bitmask_set(FlatSim(Ifilter),Summary(Ifilter).FlatNaNMap,InPar.Bit_FlatNaN,InPar.BitClass,InPar.Bit_CombFun); % Na'ama, 2018-06-06
    end
    if (~isempty(InPar.Bit_FlatLowNim))
        %FlatSim = bitmask_set(FlatSim,Summary.FlatLowNimMap,InPar.Bit_FlatLowNim,InPar.BitClass,InPar.Bit_CombFun);
        FlatSim(Ifilter) = bitmask_set(FlatSim(Ifilter),Summary(Ifilter).FlatLowNimMap,InPar.Bit_FlatLowNim,InPar.BitClass,InPar.Bit_CombFun); % Na'ama, 2018-06-06
    end
    if (~isempty(InPar.Bit_FlatHighStd))
        %FlatSim = bitmask_set(FlatSim,Summary.FlatHighStdMap,InPar.Bit_FlatHighStd,InPar.BitClass,InPar.Bit_CombFun);
        FlatSim(Ifilter) = bitmask_set(FlatSim(Ifilter),Summary(Ifilter).FlatHighStdMap,InPar.Bit_FlatHighStd,InPar.BitClass,InPar.Bit_CombFun); % Na'ama, 2018-06-06
    end
    if (~isempty(InPar.Bit_FlatLow))
        %FlatSim = bitmask_set(FlatSim,Summary.FlatLowMap,InPar.Bit_FlatLow,InPar.BitClass,InPar.Bit_CombFun);
        FlatSim(Ifilter) = bitmask_set(FlatSim(Ifilter),Summary(Ifilter).FlatLowMap,InPar.Bit_FlatLow,InPar.BitClass,InPar.Bit_CombFun); % Na'ama, 2018-06-06
    end

   
end


    
    
