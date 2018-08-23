function [Res,Count]=bitmask_find(MaskIm,BitMask,Comparison,BitType)
% Return images of logical flags indicating if bitmask was found in mask.
% Package: @MASK
% Description: Given a MASK object and a bit mask return images of logical
%              flags indicating if bitmask was found in mask.
%              The find can be done either by exact value, or by bit and.
% Input  : - MASK object.
%          - Bit mask. This can be a string of bit mask name, a cell array
%            of bit mask names, a bit mask decimal representation or a
%            vector of bit indices.
%          - Comparison method:
%            'bitand'   - Bit and operation. Default.
%            'exactval' - By exact value.
%          - Bit type. If the BitMask is numeric then this parameter is
%            used to specify if the bit mask is a vector of indices
%            ('index') or a bit mask value ('value').
%            Default is 'value'.
% Output : - A SIM array of images in which the image field contains
%            a logical image specifing which pixel satisfy the bit mask.
%          - A matrix (elemnt per MASK image) in which each value indicate
%            the number of pixels that satisfy the bit mask condition.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,Count]=bitmask_find(S(1),16,'exactval');
% Reliable: 2
%--------------------------------------------------------------------------

MaskField     = 'Mask';
%MaskDicField  = 'MaskDic';
ImageField    = 'Im';

Def.Comparison = 'bitand';  % 'bitand'|exactval' 
Def.BitType    = 'Value';   % 'value'|'index'
if (nargin==2)
    Comparison = Def.Comparison;
    BitType    = Def.BitType;
elseif (nargin==3)
    BitType    = Def.BitType;
elseif (nargin==4)
    % do nothing
else
    error('Input arguments: [MaskIm,BitMask,[Comparison,BitType]]');
end


if (ischar(BitMask))
    BitMask = {BitMask};
end


% Define output SIM array
Res   = SIM(size(MaskIm));

Nmask = numel(MaskIm);
for Imask=1:1:Nmask
    % for each mask image
    
    if (iscell(BitMask))
        % get bit map value from names
        BitMaskVal = bitname2val(MaskIm(Imask),BitMask);

    else
        switch lower(BitType)
            case 'value'
                % assume bit mask is provided by its decimal value
                BitMaskVal = BitMask;
            case 'index'
                BitMaskVal = sum(bitset(0,BitMask));
            otherwise
                error('Unown BitType option');
        end
    end

    
    
    if (isempty(MaskIm(Imask).(MaskField)))
        % Mask is empty - nothing to find
        % don't touch the Mask
    else
        switch lower(Comparison)
            case 'exactval'
                Res(Imask).(ImageField) = MaskIm(Imask).(MaskField) == BitMaskVal;
            case 'bitand'
                Res(Imask).(ImageField) = bitand(MaskIm(Imask).(MaskField),BitMaskVal)>0;
            otherwise
                error('Unknown Comparison option');
        end
    end
    
    if (nargout>1)
        Count(Imask) = sum(Res(Imask).(ImageField)(:));
    end
    
end
        
