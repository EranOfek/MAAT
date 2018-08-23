function Sim=replace_with_noise(Sim,varargin)
% Replace a range of values by noise.
% Package: @SIM
% Description: Replace some values, or ranges (e.g., 0, NaN) in SIM images
%              with the image background and noise (in pixels not in the
%              specified ranges) or by constant value.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Range' - A two column matrix of ranges [Low High].
%                      Any pixel which value within these ranges will be
%                      replaced. NaN will replace NaN values.
%                      Default is [0 0; NaN NaN].
%            'ExecField' - Cell array of SIM fields on which to execute 
%                      this function. Default is {'Im'}.
%            'Noise' - Noise type:
%                      'back'  - Replace values in ranges with the global
%                                background of the image found using
%                                mode_fit.m
%                      'const' - Replace values in ranges by a constant
%                                given by 'Const'.
%                      'norm'  - Replace values in ranges by normally
%                                distributed values which mean is the 
%                                global background and std is the global
%                                std (found using mode_fit.m).
%                      'poiss' - Replace values in ranges by Poisson
%                                distributed values which mean is the 
%                                global background
%                                (found using mode_fit.m).
%                      Default is 'norm'.
%            'Const' - A constant to use in the 'Noise'=='const' option.
%                      Default is 0.
%            'ReplaceKey' - A three column cell array of {key,val,comment},
%                           or an HEAD object which keywords to replace
%                           or add to the SIM header.
%            'AddKey'     - Like 'ReplaceKey', but adding keywords,
%                           without replacment.
%            'MaskFun'    - Function that sets the bit mask:
%                           Fun(Sim1,MaskFunPar{:}).
%                           Default is empty.
%            'MaskFunPar' - Additional parameters to pass to MaskFun.
%                           Default is {}.
% Output : - A SIM object with the replaced pixels.
% See also: SIM/replace.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=SIM; Sim.Im=[0 0 1 NaN; 10 19 15 11; 13 12 11 0];
%          Sim=replace_with_noise(Sim)
% Reliable: 2
%--------------------------------------------------------------------------


DefV.Range              = [0 0; NaN NaN];    % range of values 
DefV.ExecField          = {SIM.ImageField};
DefV.Noise              = 'norm';   % 'back'|'const'|'norm'|'poiss'
DefV.Const              = 0;
DefV.ReplaceKey         = {};
DefV.AddKey             = {};
DefV.MaskFun            = [];   % should work on SIM! and return a SIM
DefV.MaskFunPar         = {};   % e.g., MaskFun=@eq, MaskFunPar={0}
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ExecField))
    InPar.ExecField  = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

Nsim = numel(Sim);
for Isim=1:1:Nsim
    % for each SIM element
    for If=1:1:Nf
        % for each SIM field
        
        % find of ranges
        Flag = Util.array.find_ranges_flag(Sim(Isim).(InPar.ExecField{If}),InPar.Range);
        FlagIR = Flag==0;  % in range
        FlagNR = Flag>0;   % nor in range
        % Calculate the mode and std at location which were not found in
        % range
        [B,S] = Util.stat.mode_fit(Sim(Isim).(InPar.ExecField{If})(FlagIR));
        
        % Replace the pixels found in range with some value
        switch lower(InPar.Noise)
            case 'const'
                % replace by a user supplied constant number
                Sim(Isim).(InPar.ExecField{If})(FlagNR) = InPar.Const;
            case 'back'
                % replace by the constant background level
                Sim(Isim).(InPar.ExecField{If})(FlagNR) = B;
            case 'norm'
                % replace by the background with normally distributed noise
                Sim(Isim).(InPar.ExecField{If})(FlagNR) = randn(sum(FlagNR(:)),1).*S + B;
            case 'poiss'
                % replace by the background with Poisson distributed noise
                Sim(Isim).(InPar.ExecField{If})(FlagNR) = poissrnd(B,sum(FlagNR(:)),1);
            otherwise
                error('Unknown Noise option');
        end
    end
    
    %--- Update mask ---
    if (~isempty(InPar.MaskFun))
        Sim(Isim) = InPar.MaskFun(Sim(Isim),InPar.MaskFunPar{:});
    end
    
    %--- Update header ---
    % replace
    if (~isempty(InPar.ReplaceKey))
        Sim(Isim) = replace_key(Sim(Isim),InPar.ReplaceKey);
    end
    % add
    if (~isempty(InPar.AddKey))
        Sim(Isim) = add_key(Sim(Isim),InPar.AddKey);
    end
    
    
end
        
    


