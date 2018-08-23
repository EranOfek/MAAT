function MaskIm=bitmask_set(MaskIm,Flag,BitName,BitType,CombFun)
% Set on specific bits in some pixels in a MASK object.
% Package: @MASK
% Description: Set specific bits in pixels in a MASK object
%              (e.g., a SIM image).
%              This function get a MASK object, a logical image of flags
%              indicating which pixels to set, the bits to set and the
%              set function (e.g., or/and).
% Input  : - A MASK object.
%          - A logical image indicating which pixels to set (true),
%            or not to touch (false).
%          - The bits to set. This can be either a decimal value,
%            or a string or cell array of strings of bits to set
%            (bit names are defined in the dictionary).
%          - Bit type. Default is 'uint32'.
%          - Combining function: @bitor|@bitand|@bitxor. Default is @bitor.
% Output : - The input MASK object with the setted bits.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S3=bitmask_set(S2,logical(Flag),3)
% Reliable: 2
%--------------------------------------------------------------------------

ImageField   = 'Im';
MaskField    = 'Mask';

Def.BitType = 'uint32';
Def.CombFun = @bitor;
if (nargin==3)
    BitType = Def.BitType;
    CombFun = Def.CombFun;
elseif (nargin==4)
    CombFun = Def.CombFun;
elseif (nargin==5)
    % do nothing
else
    error('Illegal number of input arguments');
end


Nmask = numel(MaskIm);
if (SIM.issim(Flag))
    FlagIm = Flag;
else
    FlagIm = SIM;
    FlagIm.(ImageField) = Flag;
end
Nflag = numel(FlagIm);

if (~isnumeric(BitName))
    Value = bitname2val(MaskIm,BitName);
else
    Value = BitName;
end

if (isnan(Value))
    error('Unknown Bit name');
end

N = max(Nmask,Nflag);
for I=1:1:N
    % for each image
    Imask = min(I,Nmask);
    Iflag = min(I,Nflag);
    
    if (isempty(MaskIm(Imask).(MaskField)))
        MaskIm(Imask).(MaskField) = zeros(size(FlagIm(Iflag).(ImageField)),BitType);
    end
        
    MaskIm(Imask).(MaskField)(FlagIm(Iflag).(ImageField)) = CombFun(MaskIm(Imask).(MaskField)(FlagIm(Iflag).(ImageField)),Value(Iflag));
end


    



