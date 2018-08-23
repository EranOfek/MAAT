function DicMap=showdic(MaskIm,DefDic)
% Show MASK object dictionary.
% Package: @MASK
% Description: Get dictionary for a MASK object and print it to screen.
% Input  : - MASK object.
%          - Default dictionary function name or structure containing
%            dictionary. Defauly is @MASK.def_bitmask_pipeline.
% Output : - Dictionary map.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: DicMap=showdic(S);
%          [~,~,BitSt]=val2bitname(S,17);showdic(BitSt)
% Reliable: 2
%--------------------------------------------------------------------------

MaskDicField  = 'MaskDic';

Def.Dic       = @MASK.def_bitmask_pipeline; 
if (nargin==1),
    DefDic = Def.Dic;
elseif (nargin==2),
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isempty(MaskIm(1).(MaskDicField))),
    Dic = DefDic;
else
    Dic = MaskIm(1).(MaskDicField);
end

if (isa(Dic,'function_handle')),
    [~,DicMap] = Dic();
else
    DicMap = Dic;
end

Ndic = numel(DicMap);
fprintf('  Bit Name                Index    Decimal\n')
fprintf('  ----------------------  -----    -------\n')
for Idic=1:1:Ndic,
    fprintf('  %-25s  %2d    %d\n',DicMap(Idic).Name,DicMap(Idic).Ind,2.^DicMap(Idic).Ind);
end

    