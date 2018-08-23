function [Sim,SimFlag]=flag_backgrad(Sim,varargin)
% Flag high values of SIM images background gradient and add to MASK
% Package: @SIM
% Description: Given a SIM object, for each image background flag values
%              with large background gradient in the SIM bit mask.
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


BackField    = SIM.BackField;
ErrField     = SIM.ErrField;

DefV.Nsigma              = 5;
DefV.Bit_BackGrad        = 'Bit_BackGrad';
DefV.MaskDic             = @MASK.def_bitmask_pipeline;

DefV.BackPar             = {};


DefV.ExecField           = SIM.ImageField;
DefV.BitType             = 'uint32';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (nargout>1)
    SimFlag = SIM(size(Sim));
end

% check if background exist - if not then generate background
IsBackPop = isfield_populated(Sim,'BackField') & isfield_populated(Sim,'ErrField');
Sim(~IsBackPop) = background(Sim(~IsBackPop),InPar.BackPar{:});


Nsim = numel(Sim);
for Isim=1:1:Nsim
    [GX,GY] = gradient(Sim(Isim).(BackField));
    Grad = sqrt(GX.^2 + GY.^2);
    Grad./Sim(Isim).(ErrField);
    
    Tmp = (Sim(Isim).(InPar.ExecField) - Sim(Isim).(BackField))./Sim(Isim).(ErrField) < -InPar.Nsigma;

    Sim(Isim) = bitmask_set(Sim(Isim),Tmp,InPar.Bit_Hole,InPar.BitType);

    if (nargout>1)
        Sim(Isim),SimFlag(Isim).(InParExecField) = Tmp;
    end
end