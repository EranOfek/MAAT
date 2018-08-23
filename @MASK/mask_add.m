function MaskIm1=mask_add(MaskIm1,MaskIm2,CombFun)
% Combine the corresponding Mask images of two MASK objects.
% Package: @MASK
% Description: Combine the corresponding Mask images of two MASK
%              object arrays and store them in the first object array.
% Input  : - First MASK object.
%          - Second MASK object. Number of elements should be identical to
%            that in the first input object or 1.
%          - Combine function: @bitand, @bitor. Default is @bitor.
% Output : - The first MASK object, in which the Mask field is combined
%            with that of the second MASK object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S3=mask_add(S1,S2);
% Reliable: 2
%--------------------------------------------------------------------------

    
MaskField      = 'Mask';
%MaskDicField   = 'MaskDic';

if (nargin<3)
    CombFun = @bitor;
end
FunName = func2str(CombFun);
    
Nm1 = numel(MaskIm1);
Nm2 = numel(MaskIm2);
if ~(Nm1==Nm2 || Nm2==1 || Nm2==0)
    error('Illegal number of elements in the second Mask image');
end

N   = max(Nm1,Nm2);

for I=1:1:N
    Im1 = min(I,Nm1);
    Im2 = min(I,Nm2);
    
    if (isempty(MaskIm1(Im1).(MaskField)) && ~isempty(MaskIm2(Im2).(MaskField))),
        % Mask1 is empty, Mask2 is not empty
        switch lower(FunName)
            case 'bitor'
                MaskIm1(Im1).(MaskField) = MaskIm2(Im2).(MaskField);
            case 'bitand'
                % do nothing
            case 'bitxor'
                MaskIm1(Im1).(MaskField) = zeros(size(MaskIm2(Im2).(MaskField)));
                MaskIm1(Im1).(MaskField) = CombFun(MaskIm1(Im1).(MaskField),MaskIm2(Im2).(MaskField));
            otherwise
                error('Unknown FunName option');
        end
    elseif (~isempty(MaskIm1(Im1).(MaskField)) && isempty(MaskIm2(Im2).(MaskField))),
        % Mask1 is not empty, Mask2 is empty      
        switch lower(FunName)
            case 'bitor'
                MaskIm1(Im1).(MaskField) = MaskIm1(Im1).(MaskField);
            case 'bitand'
                % do nothing
            case 'bitxor'
                MaskIm2(Im2).(MaskField) = zeros(size(MaskIm1(Im1).(MaskField)));
                MaskIm1(Im1).(MaskField) = CombFun(MaskIm1(Im1).(MaskField),MaskIm2(Im2).(MaskField));
            otherwise
                error('Unknown FunName option');
        end  
    else
        MaskIm1(Im1).(MaskField) = CombFun(MaskIm1(Im1).(MaskField),MaskIm2(Im2).(MaskField));
    end
    
end
    
