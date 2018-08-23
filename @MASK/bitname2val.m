function [Value,BitInd]=bitname2val(MaskIm,BitName,Dic,DefMaskClass)
% List of bitmask names to a decimal bitmask value.
% Package: @MASK
% Description: Convert a list of bitmask names to a decimal bitmask value,
%              and bit indices.
% Input  : - MASK object.
%          - String of bitmask name or a cell array of bitmask names.
%          - Bitmap dictionary. By default, the MASK object dictionary will
%            be used. If the MASK object does not contain a dictionary
%            the use the default dictionary (@MASK.def_bitmask_pipeline).
%            Dictionary can be either a function handle or a structure
%            array.
%          - Mask class. By default will use the class of the MASK image.
%            If mask image does not exits then default is 'uint32'.
% Output : - A vector (element per mask image) of bitmask values. Each
%            value is the decimal representation of all the bitmask names.
%          - A matrix of bit indices. Column per mask image, row per bit.
%            0 if bit name was not in dictionary.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:  Value=bitname2val(S,'Bit_Edge')
%           Value=bitname2val(S,{'Bit_Edge','Bit_ImSaturated'})
% Reliable: 2
%--------------------------------------------------------------------------


MaskField     = 'Mask';
MaskDicField  = 'MaskDic';


Def.Dic       = @MASK.def_bitmask_pipeline; 
Def.MaskClass = 'uint32';
if (nargin==2)
    Dic          = Def.Dic;
    DefMaskClass = Def.MaskClass;
elseif (nargin==3)
    DefMaskClass = Def.MaskClass;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end


if (~iscell(BitName))
    BitName = {BitName};
end
Nbitname = numel(BitName);


Nmask = numel(MaskIm);
Value = zeros(size(MaskIm));
BitInd = zeros(Nbitname,Nmask);
for Imask=1:1:Nmask
    if (isempty(MaskIm(Imask).(MaskDicField)))
        % Mask dictonary is not available in Mask image
        % use default dicitonary
        DicI = Dic;
    else
        DicI = MaskIm(Imask).(MaskDicField);
    end
    
    if (isempty(MaskIm(Imask).(MaskField)))
        MaskClass = DefMaskClass;
    else
        MaskClass = class(MaskIm(Imask).(MaskField));
    end
    
    % get the dictionary map
    if (isa(DicI,'function_handle'))
        % Dictionary is a function
        [~,DicMap] = DicI();
    else
        % assume dictionary is a structure
        DicMap = DicI;
    end
    
    %DicInd = [DicMap.Ind];
    Ind = zeros(Nbitname,1);
    for Ibitname=1:1:Nbitname
        IndName = strcmpi({DicMap.Name},BitName{Ibitname});
        
        if (isempty(Ind))
            warning('BitName %s was not found in dictionary',BitName{Ibitname});
        else
            if (sum(IndName)~=0)
                Ind(Ibitname) = find(IndName);
            end
        end
    end
    BitInd(:,Imask) = Ind;
    
    Ind = Ind(Ind>0);
    Value(Imask) = sum(bitset(zeros(1,1,MaskClass),Ind));
end
    
% replace 0 by NaN
Value(Value==0) = NaN;

