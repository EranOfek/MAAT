function Sim=interp2_overpix(Sim,varargin)
%--------------------------------------------------------------------------
% interp2_overpix function                                      class/@SIM
% Description: Interp images in SIM array over bad pixels.
%              The bad pixels are specified by flag image or specific
%              values in the bit mask.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField'  - Fields in the SIM object over which to
%                           interpolate bad pixels.
%                           Default is {'Im'}.
%            'Method'     - inpaint_nans.m interpolation method.
%                           Default is 0.
%            'FlagMap'    - Map of logical flags. Pixels which value is
%                           true will be interpolated over.
%                           This can be either a matrix of logicals,
%                           a cell array of logicals, or a SIM object
%                           in which the 'Im' field is logicals.
%                           If emptym then will use the 'Mask' field
%                           in the SIM object to define the bad pixels.
%                           Default is empty.
%            'Bit_Over2Interp'- A numeric value of bit mask, or a cell
%                           array of bit mask names to interpolate.
%                           Pixels for which this specific bits in the
%                           'Mask' fields are on, will be interpolated
%                           over.
%                           Default is
%                           {'Bit_BadPixel','Bit_BiasNoise0','Bit_FlatNaN','Bit_CR'}.
%            'MaskDic'    - Default mask dictionary
%            'Bit_PixInterp'- The bit mask value or name of the bit in the
%                           bit mask which to set to on for interpolated
%                           pixels. Default is 'Bit_PixInterp'.
%            'BitClass'   - Default bit mask class. Default is 'uint32'.
% Output : - A SIM object with the interpolated pixels and updated mask.
% See also: inpaint_nans.m, SIM/interp2_nan.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=interp2_overpix(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField        = 'Im';
MaskDicField      = 'MaskDic';

DefV.ExecField          = {ImageField};
DefV.Method             = 0;
DefV.FlagMap            = [];
DefV.Bit_Over2Interp    = {'Bit_BadPixel','Bit_BiasNoise0','Bit_FlatNaN','Bit_CR'};
DefV.MaskDic            = @MASK.def_bitmask_pipeline;
DefV.Bit_PixInterp      = 'Bit_PixInterp';
DefV.BitClass           = 'uint32';

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end

Nf   = numel(InPar.ExecField);
Nsim = numel(Sim);
for Isim=1:1:Nsim
    if (isempty(InPar.FlagMap))
        % interpolate over pixels defined by bit mask
        % look for pixels
        % check that MaskDic exist
        if (isempty(Sim(Isim).(MaskDicField)))
            Sim(Isim).(MaskDicField) = InPar.MaskDic;
        end
        
        % generate a matrix of logical flags of pixels over which to
        % interpolate the images
        Res = bitmask_find(Sim(Isim),InPar.Bit_Over2Interp,'bitand','value');
        
        FlagMap = Res.(ImageField);
    else
        % pixels to interpolate are given in FlagMap
        if (islogical(InPar.FlagMap))
            FlagMap = InPar.FlagMap;
        elseif (iscell(InPar.FlagMap))
            FlagMap = InPar.FlagMap{Isim};
        elseif (SIM.issim(InPar.FlagMap))
            FlagMap = Res.(ImageField);
        else
            error('Unknown FlagMap option');
        end
    end
    
    for If=1:1:Nf
        % set pix value to NaN
        Sim(Isim).(InPar.ExecField{If})(FlagMap) = NaN;
        % interpolate over NaN
        Sim(Isim).(InPar.ExecField{If}) = inpaint_nans(Sim(Isim).(InPar.ExecField{If}),InPar.Method);
    end
    
    % Set the interpolation bit mask
    Sim(Isim) = bitmask_set(Sim(Isim),FlagMap,InPar.Bit_PixInterp,InPar.BitClass,@bitor);
    
end


        