function [BitName,BinInd,BitSt]=val2bitname(MaskIm,Value,Dic,varargin)
% Convert a bit mask value to a list of flagged bit names.
% Package: @MASK
% Description: Convert a bit mask value to a list of flagged bit names.
% Input  : - MASK object.
%          - Decimal representation of the bitmask.
%          - Bitmap dictionary. By default, the MASK object dictionary will
%            be used. If the MASK object does not contain a dictionary
%            the use the default dictionary (@MASK.def_bitmask_pipeline).
%            Dictionary can be either a function handle or a structure
%            array.
% Output : - Cell array of bit names.
%          - Vector of bit indices.
%          - Structure array of bit names and indices.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [BitName,BinInd]=val2bitname(S,16);
%          [BitName,BinInd]=val2bitname(S,17);
% Reliable: 2
%--------------------------------------------------------------------------


%MaskField     = 'Mask';
MaskDicField  = 'MaskDic';


Def.Dic       = @MASK.def_bitmask_pipeline; 
if (nargin==2)
    Dic          = Def.Dic;
elseif (nargin==3)
    % do nothing
else
    error('Illegal number of input arguments');
end

BinInd = strfind(fliplr(dec2bin(Value)),'1');

Nmask = numel(MaskIm);
if (Nmask>1)
    warning('Return bitname only for first Mask image');
end

    Imask = 1;    
    if (isempty(MaskIm(Imask).(MaskDicField)))
        % Mask dictonary is not available in Mask image
        % use default dicitonary
        DicI = Dic;
    else
        DicI = MaskIm(Imask).(MaskDicField);
    end
    
     % get the dictionary map
    if (isa(DicI,'function_handle'))
        % Dictionary is a function
        [~,DicMap] = DicI();
    else
        % assume dictionary is a structure
        DicMap = DicI;
    end
    
    DicInd = Util.array.findmany([DicMap.Ind],BinInd);
    
    BitSt = DicMap(DicInd);
    BitName = {BitSt.Name};
    
