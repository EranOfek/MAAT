function Sim=bfun2sim(Sim1,Sim2,Op,varargin)
% Run a bi-parameter function (e.g., @plus) on two SIM arrays.
% Package: @SIM
% Description: Run a bi-parameter function (e.g., @plus) on two SIM
%              image objects, or a SIM object and another input.
% Input  : - The first SIM object.
%          - The second parameter. This can be one of the following:
%            A SIM object,; a cell array of images; an array of scalars;
%            a function handle that returns an image (matrix);
%            a PSF object, or empty. If empty then will attempt to read
%            the PSF from the first SIM input.
%            The number of elements is either 1 or equal to the number
%            of elements in the first SIM input.
%          - Function handle for the operator (e.g., @plus).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FunAddPar' - A cell array of additional parameters to pass
%                          to the operator function. Default is {}.
%            'Sim2FunPar'- If the second input argument is a function
%                          handle that is used to calculate the second
%                          image, then this is a cell array of arguments
%                          to pass to this function. Default is {}.
%            'ExecField' - The SIM fields on which to execute the operator.
%                          Default is {'Im'}.
%            'CCDSEC1'   - CCDSEC for the first input argument.
%                          This can be a string of the header keyword,
%                          or [Xmin, Xmax, Ymin, Ymax].
%                          If empty use entire image. Default is [].
%            'CCDSEC2'   - CCDSEC for the second input argument.
%                          Default is [].
%            'HeadFrom'  - From which input argument to copy the header
%                          from: 1|2|[]. If empty do not copy header.
%                          Default is 1.
%            'CatFrom'   - From which input argument to copy the Cat field
%                          from: 1|2|[]. If empty do not copy Cat.
%                          Default is 1.
%            'ReplaceKey'- A three column cell array of {key,val,comment},
%                          or an HEAD object which keywords to replace
%                          or add to the SIM header.
%                          Default is {}.
%            'AddKey'    - Like 'ReplaceKey', but adding keywords,
%                          without replacment.
%                          Default is {}.
%            'MaskFun'    - Function that sets the bit mask:
%                           Fun(Sim1,Sim2,MaskFunPar{:}).
%                           Default is empty.
%            'MaskFunPar' - Additional parameters to pass to MaskFun.
%                           Default is {}.
% Output : - A SIM object with the output operation.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S1=bfun2sim(S1,S2,@plus);
%          S1=bfun2sim(S1,3,@minus);
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField  = HEAD.HeaderField;
CatField     = AstCat.CatField;

DefV.FunAddPar            = {};
DefV.Sim2FunPar           = {};
DefV.ExecField            = {SIM.ImageField};
DefV.CCDSEC1              = [];
DefV.CCDSEC2              = [];
DefV.HeadFrom             = 1;   % 1|2|[]
DefV.CatFrom              = 1;  % 1|2|[]
DefV.ReplaceKey           = {};
DefV.AddKey               = {};
DefV.MaskFun              = [];   % should work on SIM! and return a SIM
DefV.MaskFunPar           = {};   % e.g., MaskFun=@eq, MaskFunPar={0}

if (isempty(varargin))
    InPar = DefV;
else
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
    %InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
end


if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

Nsim1 = numel(Sim1);
Nsim2 = numel(Sim2);
if ~(Nsim1==Nsim2 || Nsim1==1 || Nsim2==1)
    error('The number of elements in the input SIM objects should be equal or 1');
end
N     = max(Nsim1,Nsim2);

Sim = SIM(N,1);
for I=1:1:N
    % for each image
    I1 = min(I,Nsim1);
    I2 = min(I,Nsim2);
    
    for If=1:1:Nf
        if (~isempty(Sim1(I1).(InPar.ExecField{If})))
            % ExecField in Sim1 is not empty:
            if (SIM.issim(Sim2))
                Tmp2 = Sim2(I2).(InPar.ExecField{If});
            elseif (iscell(Sim2))
                Tmp2 = Sim2{I2};
            elseif (isnumeric(Sim2))
                % Input is vector of scalars
                % in order to provide a matrix use the cell input...
                Tmp2 = Sim2(I2);
            elseif (isa(Sim2,'function_handle'))
                % Sim2 is a function
                Tmp2 = Sim2(InPar.Sim2FunPar{:});
            elseif (ClassPSF.isClassPSF(Sim2))
                Tmp2 = getpsf(Sim2(I2));
            elseif (isempty(Sim2))
                error('Empty Sim2 not yet available');
                % Attempt to read the PSF field from the first SIM input
            else
                error('Unknown Sim2 option');
            end
            if (~isempty(InPar.CCDSEC2))
                CCDSEC2 = ccdsec(Sim2(I2),InPar.CCDSEC2);
                Tmp = Tmp(CCDSEC2(3):CCDSEC2(4),CCDSEC2(1):CCDSEC2(2));
            end
            
            % Match the precision of Sim1 and Sim2 (take the highest)
            Tmp1 = Sim1(I1).(InPar.ExecField{If});
            if isinteger(Tmp1) && isfloat(Tmp2)
                Tmp1 = cast(Tmp1, 'like', Tmp2);
            elseif isfloat(Tmp1) && isinteger(Tmp2)
                Tmp2 = cast(Tmp2, 'like', Tmp1);
            end

            % ccdsec
            if (isempty(InPar.CCDSEC1))
                % operation
                Sim(I).(InPar.ExecField{If}) = Op(Tmp1,Tmp2,InPar.FunAddPar{:});
            else
                % operation with CCDSEC
                CCDSEC1 = ccdsec(Sim1(I1),InPar.CCSEC1);
                Sim(I).(InPar.ExecField{If}) = Op(Tmp1(CCDSEC1(3):CCDSEC1(4),CCDSEC1(1):CCDSEC1(2)),Tmp2,InPar.FunAddPar{:});
            end

        end
    end
    
    % catalog
    if (isempty(InPar.CatFrom))
        Sim(I).(CatField) = [];
    else
        if (InPar.CatFrom==2)
            Sim(I).(CatField) = Sim2(I2).(CatField);
        end
    end

    % header
    if (isempty(InPar.HeadFrom))
        Sim(I).(HeaderField) = [];
    else
        if (InPar.HeadFrom==2)
            Sim(I).(HeaderField) = Sim2(I2).(HeaderField);
            Sim(I).WCS = Sim2(I2).WCS; % Na'ama, 20180525
        end
        if (InPar.HeadFrom==1)
            Sim(I).(HeaderField) = Sim1(I2).(HeaderField);
            Sim(I).WCS = Sim1(I2).WCS; % Na'ama, 20180525
        end
    end

    %--- copy mask ---
    Sim(I).Mask    = Sim1(I1).Mask;
    Sim(I).MaskDic = Sim1(I1).MaskDic;
    
    %--- Update mask ---
    if (~isempty(InPar.MaskFun))
        %Sim(I) = InPar.MaskFun(Sim1(I1),Sim2(I2),InPar.MaskFunPar{:});
        Sim(I) = InPar.MaskFun(Sim(I),Sim2(I2),InPar.MaskFunPar{:}); % Na'ama, 20171205
    end
    
    %--- Update header ---
    % replace
    if (~isempty(InPar.ReplaceKey))
        %Sim(I) = replace_key(Sim1(I1),InPar.ReplaceKey);
        Sim(I) = replace_key(Sim(I),InPar.ReplaceKey); % Na'ama, 20171205
    end
    % add
    if (~isempty(InPar.AddKey))
        %Sim(I) = add_key(Sim1(I1),InPar.AddKey);
        Sim(I) = add_key(Sim(I),InPar.AddKey); % Na'ama, 20171205
    end
    
end
