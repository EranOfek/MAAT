function Res=mask_combine(MaskIm,CombFun)
% Combine all the Mask images in a single MASK object array.
% Package: @MASK
% Description: Combine all the Mask images in a single MASK object array.   
% Input  : - MASK object.
%          - Combine function: @bitand, @bitor. Default is @bitor.
% Output : - A Mask object with a single elemnt containing the combination
%            of all the elements in the input MASK object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=mask_combine(MaskIm)
% Reliable: 2
%--------------------------------------------------------------------------

MaskField      = 'Mask';
MaskDicField   = 'MaskDic';

if (nargin<2),
    CombFun = @bitor;
end

    
Nmask = numel(MaskIm);

Res = MASK;
for Imask=1:1:Nmask,
    if (isempty(Res.(MaskField))),
        % Res was empty -
        Res.(MaskField)    = MaskIm(Imask).(MaskField);
        Res.(MaskDicField) = MaskIm(Imask).(MaskDicField);
    else
        if (~isempty(MaskIm(Imask).(MaskField))),
            Res.(MaskField) = CombFun(Res.(MaskField),MaskIm(Imask).(MaskField));
        end
    end
 end


    
