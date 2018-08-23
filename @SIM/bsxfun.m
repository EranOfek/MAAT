function Sim1=bsxfun(Sim1,Sim2,Op,varargin)
% bsxfun function for SIM objects. 
% Package: @SIM
% Description: bsxfun function for SIM objects. 
% Input  : - The first SIM object.
%          - The second parameter. This can be one of the following:
%            A SIM object; a cell array of images; A vector.
%            For SIM object and cell array input the number of elements
%            is either 1 or equal to the number of elements in the first
%            SIM input. Each element (or the input vector) is a vector
%          - Function handle for the operator (e.g., @plus).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - The SIM fields on which to execute the operator.
%                          Default is {'Im'}.
%            'CCDSEC1'   - CCDSEC for the first input argument.
%                          This can be a string of the header keyword,
%                          or [Xmin, Xmax, Ymin, Ymax].
%                          If empty use entire image. Default is [].
%            'HeadFrom'  - From which input argument to copy the header
%                          from: 1|2|[]. If empty do not copy header.
%                          Default is 1.
%            'CatFrom'   - From which input argument to copy the Cat field
%                          from: 1|2|[]. If empty do not copy Cat.
%                          Default is [].
%            'ReplaceKy'- A three column cell array of {key,val,comment},
%                          or an HEAD object which keywords to replace
%                          or add to the SIM header.
%            'AddKey'    - Like 'ReplaceKey', but adding keywords,
%                          without replacment.
%            'MaskFun'    - Function that sets the bit mask:
%                           Fun(Sim1,MaskFunPar{:}).
%                           Default is empty.
%            'MaskFunPar' - Additional parameters to pass to MaskFun.
%                           Default is {}.
% Output : - A SIM object with the output operation.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=bsxfun(S(1),ones(4096,1),@plus);
% Reliable: 2
%--------------------------------------------------------------------------


ImageField   = 'Im';
HeaderField  = 'Header';
CatField     = 'Cat';


%DefV.FunAddPar            = {};
DefV.ExecField            = {ImageField};
DefV.CCDSEC1              = [];
%DefV.CCDSEC2              = [];
DefV.HeadFrom             = 1;   % 1|2|[]
DefV.CatFrom              = [];  % 1|2|[]
DefV.ReplaceKey           = {};
DefV.AddKey               = {};
DefV.MaskFun              = [];   % should work on SIM! and return a SIM
DefV.MaskFunPar           = {};   % e.g., MaskFun=@eq, MaskFunPar={0}
if (isempty(varargin))
    InPar = DefV;
else
    %InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
end


if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

N1 = numel(Sim1);
N2 = numel(Sim2);

for I1=1:1:N1
    % for each image
    
    for If=1:1:Nf
        if (SIM.issim(Sim2))
            I2  = min(I1,N2);
            Tmp = Sim2(I2).(InPar.ExecField{If});
        elseif (iscell(Sim2))
            I2  = min(I1,N2);
            Tmp = Sim2{I2};
        else
            Tmp = Sim2;
        end
        % no ccdsec for Sim2
        
        % ccdsec for Sim1
        if (isempty(InPar.CCDSEC1))
            % operation
            Sim1(I1).(InPar.ExecField{If}) = bsxfun(Op,Sim1(I1).(InPar.ExecField{If}),Tmp);
        else
            % operation with CCDSEC
            CCDSEC1 = ccdsec(Sim1(I1),InPar.CCSEC1);
            Sim1(I1).(InPar.ExecField{If}) = bsxfun(Op,Sim1(I1).(InPar.ExecField{If})(CCDSEC1(3):CCDSEC1(4),CCDSEC1(1):CCDSEC1(2)),Tmp);
        end
    end
    
    % catalog
    if (isempty(InPar.CatFrom))
        Sim1(I1).(CatField) = [];
    else
        if (InPar.CatFrom==2)
            Sim1(I1).(CatField) = Sim2(I2).(CatField);
        end
    end

    % header
    if (isempty(InPar.HeadFrom))
        Sim1(I1).(HeaderField) = [];
    else
        if (InPar.HeadFrom==2)
            Sim1(I1).(HeaderField) = Sim2(I2).(HeaderField);
        end
    end

    %--- Update mask ---
    if (~isempty(InPar.MaskFun))
        Sim1(I2) = InPar.MaskFun(Sim1(I1),InPar.MaskFunPar{:});
    end
    
    %--- Update header ---
    % replace
    if (~isempty(InPar.ReplaceKey))
        Sim1(I1) = replace_key(Sim1(I1),InPar.ReplaceKey);
    end
    % add
    if (~isempty(InPar.AddKey))
        Sim1(I1) = add_key(Sim1(I1),InPar.AddKey);
    end
    
end
