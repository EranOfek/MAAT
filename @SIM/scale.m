function [Sim,Scale]=scale(Sim,Scale,varargin)
% Scale the images in a SIM obect.
% Package: @SIM
% Description: Scale the images in a SIM obect array by a scalar
%              multiplication or any other binary function.
%              The scaling can be a vector of constant,
%              or some functions on the images.
% Input  : - A SIM object array.
%          - Scale factor. This can be one of the following:
%            A scalar that will multiply all images.
%            A vector that each element will multiply the corresponding
%            image in the SIM array.
%            A cell array of scalars (element per image).
%            A function handle: Fun(Sim,FunAddPar{:}), that gets a SIM
%            array and return a scalar per SIM image (e.g., using
%            ufun2scalar.m).
%            A string containing one of the following options calculated
%            from the 'Im' field of the SIM array:
%            'none'   - do nothing. Scale is 1.
%            'mean'   - scale each image by its nanmean
%            '1/mean' - scale each image by its 1/nanmean.
%            'median' - scale each image by its nanmedian.
%            '1/median'-scale each image by its 1/nanmedian.
%            'min'    - scale each image by its min.
%            '1/min'  - scale each image by its 1/min.
%            'max'    - scale each image by its max.
%            '1/max'  - scale each image by its 1/max.
%            'sum'    - scale each image by its sum.
%            '1/sum'  - scale each image by its 1/sum.
%            'gain'   - scale each image by its gain obtained using
%                       getkey_gain.m.
%            '1/gain' - scale each image by its 1/gain obtained using
%                       getkey_gain.m.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ScaleFun'  - The function for scaling. Default is @times.
%                          E.g., use @minus in order to subtract a constant
%                          from each image.
%            'ExecField' - The SIM fields on which to execute the operator.
%                          Default is {'Im'}.
%            'CCDSEC'    - CCDSEC for the input argument.
%                          This can be a string of the header keyword,
%                          or [Xmin, Xmax, Ymin, Ymax].
%                          If empty use entire image. Default is [].
%            'FunAddPar' - A cell array of arguments to pass to the
%                          scale function input if it is a function handle.
%                          Default is {}.
%            'GetGainPar'- A cell array of arguments to pass to the
%                          getkey_gain.m function. Default is {}.
% Output : - A SIM object array with the scaled images.
%          - An array of the scaling (multiplication) factors applied
%            to each image.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: scale(Sim,'@maen);
%          [Sim,Sc] = scale(Sim,'median');
%          [Sim,Sc] = scale(Sim,'1/gain');
%          [Sim,Sc] = scale(Sim,rand(size(Sim)));
% Reliable: 2
%--------------------------------------------------------------------------

DefV.ScaleFun           = @times;
DefV.ExecField          = {'Im'};
DefV.CCDSEC             = [];
DefV.FunAddPar          = {};
DefV.GetGainPar         = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end

Nsim = numel(Sim);
Nf   = numel(InPar.ExecField);

if (isnumeric(Scale))
    if (numel(Scale)==1)
        Scale = Scale.*ones(size(Sim));
    end
elseif (iscell(Scale))
    Scale = cell2mat(Scale);
elseif (ischar(Scale))
    switch lower(Scale)
        case 'none'
            Scale = ones(size(Sim));
        case 'mean'
            Scale = nanmean(Sim);
        case '1/mean'
            Scale = 1./nanmean(Sim);    
        case 'median'
            Scale = nanmedian(Sim);
        case '1/median'
            Scale = 1./nanmedian(Sim);
        case 'min'
            Scale = min(Sim);
        case '1/min'
            Scale = 1./min(Sim);
        case 'max'
            Scale = max(Sim);
        case '1/max'
            Scale = 1./max(Sim);
        case 'sum'
            Scale = sum(Sim);
        case '1/sum'
            Scale = 1./sum(Sim);    
        case 'gain'
            % get GAIN from header
            Scale = getkey_gain(Sim,InPar.GetGainPar{:});
        case '1/gain'
            Scale = 1./getkey_gain(Sim,InPar.GetGainPar{:});
        otherwise
            error('Unknown Scale method option');
    end
elseif (isa(Scale,'function_handle'))
    Scale = Scale(Sim,InPar.FunAddPar{:});
else
    error('Unknwon Scale option');
end

   

if (any(Scale~=1))
    for Isim=1:1:Nsim
        % for each SIM element
        for If=1:1:Nf
            % for each SIM element field

            if (isempty(InPar.CCDSEC))
%                 Sim(Isim).(InPar.ExecField{If}) = InPar.ScaleFun(Sim(Isim).(InPar.ExecField{If}),Scale(If));
                switch char(InPar.ScaleFun) % Na'ama, 20180523
                    case {'times', 'ldivide', 'rdivide'}
                        Sim(Isim).(InPar.ExecField{If}) = InPar.ScaleFun(double(Sim(Isim).(InPar.ExecField{If})),Scale(Isim));
                    otherwise
                        Sim(Isim).(InPar.ExecField{If}) = InPar.ScaleFun(Sim(Isim).(InPar.ExecField{If}),Scale(Isim));
                end
            else
                CCDSEC = ccdsec(Sim(Isim),InPar.CCDSEC);
%                 Sim(Isim).(InPar.ExecField{If})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2)) = ...
%                                     InPar.ScaleFun(Sim(Isim).(InPar.ExecField{If})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2)),Scale(If));
                switch char(InPar.ScaleFun) % Na'ama, 20180523
                    case {'times', 'ldivide', 'rdivide'}
                        Sim(Isim).(InPar.ExecField{If})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2)) = ...
                                    InPar.ScaleFun(double(Sim(Isim).(InPar.ExecField{If})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2))),Scale(Isim));
                    otherwise
                        Sim(Isim).(InPar.ExecField{If})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2)) = ...
                                    InPar.ScaleFun(Sim(Isim).(InPar.ExecField{If})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2)),Scale(Isim));
                end
            end
            
            % old
%             if (isempty(InPar.CCDSEC))
%                 Sim(Isim).(InPar.ExecField{If}) = InPar.ScaleFun(Sim(Isim).(InPar.ExecField{If}),Scale(If));
%             else
%                 CCDSEC = ccdsec(Sim(Isim),InPar.CCDSEC);
%                 Sim(Isim).(InPar.ExecField{If})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2)) = ...
%                                     InPar.ScaleFun(Sim(Isim).(InPar.ExecField{If})(CCDSEC(3):CCDSEC(4),CCDSEC(1):CCDSEC(2)),Scale(If));
%             end
        end
    end   
end
