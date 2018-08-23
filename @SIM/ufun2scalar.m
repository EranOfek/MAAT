function varargout=ufun2scalar(Sim,Fun,varargin)
% Run a scalar output uni-parameter function on a SIM object.
% Package: @SIM
% Description: Run a scalar output uni-parameter function (e.g., meannd)
%              on a SIM object, and return the output in a matrix
%              of scalar per image.
% Input  : - A SIM object
%          - A function handle (e.g., @fft).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FunAddPar'  - A cell array of additional parameters to pass
%                           to the function. Default is {}.
%            'ExecField'  - A cell array (or a single string) of field
%                           names on which to execute the function.
%                           Default is {'Im'}.
%            'CCDSEC'     - CCD section on which to  execute the function.
%                           This is either a string indicating the header
%                           keyword containing the CCDSEC or a vector of
%                           [Xmin, Xmax, Ymin, Ymax].
%            'RunOnVec'   - Run the function on the SIM image as a vector
%                           (e.g., mean(Image(:), contrary to mean(Image)).
%                           {true|false}. Default is true.
% Output : * Arbitrary number of output arguments. Each output argument
%            corresponds to a field in the 'ExecField' list (i.e., the
%            default is the single 'Im' field).
%            Each element in each output argument corresponds to one
%            elelemt (i.e., an image) in the SIM array.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: M=ufun2scalar(Sim,@mean);
%          M=ufun2scalar(Sim,@meannd,'RunOnVec',false);
% Reliable: 2
%--------------------------------------------------------------------------


ImageField      = 'Im';

DefV.FunAddPar            = {};
DefV.ExecField            = {ImageField};
DefV.CCDSEC               = [];
DefV.RunOnVec             = true;
if (isempty(varargin))
    InPar = DefV;
else
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
end

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end

Nfn  = numel(InPar.ExecField);
Nsim = numel(Sim);
% initialize
Narg = max(Nfn,nargout);
varargout = cell(Narg,1);
for Ifn=1:1:Narg
    varargout{Ifn} = zeros(size(Sim));
end

for Isim=1:1:Nsim
    % for each SIM element
    %--- ccdsec ---
    if (isempty(InPar.CCDSEC))
        % need to expedite run time
        CCDSEC = [];
    else
        CCDSEC = ccdsec(Sim(Isim));
    end
    
    %--- Uni-parameter function ---
    for Ifn=1:1:Narg
        % for each field 
        if (InPar.RunOnVec)
            % run function on SIM image as a vector
            if (isempty(CCDSEC))
                varargout{Ifn}(Isim) = Fun(Sim(Isim).(InPar.ExecField{Ifn})(:),InPar.FunAddPar{:});
            else
                Tmp = Sim(Isim).(InPar.ExecField{Ifn})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2));
                varargout{Ifn}(Isim) = Fun(Tmp(:),InPar.FunAddPar{:});
            end
        else
            % run function on SIM image as is
            if (isempty(CCDSEC))
                varargout{Ifn}(Isim) = Fun(Sim(Isim).(InPar.ExecField{Ifn}),InPar.FunAddPar{:});
            else
                varargout{Ifn}(Isim) = Fun(Sim(Isim).(InPar.ExecField{Ifn})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2)),InPar.FunAddPar{:});
            end
        end
    end
    
end

