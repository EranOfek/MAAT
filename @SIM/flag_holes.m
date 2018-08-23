function [Sim,SimFlag]=flag_holes(Sim,varargin)
% Flag holes (pixels well below background level) in SIM images
% Package: @SIM
% Description: Flag holes (pixels well below background level) in SIM
%              images.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Nsigma'  - Number of sigma below background noise for hole
%                        detection. Default is 5.
%            'Bit_Hole'- Bit hole name. Default is 'Bit_Hole'.
%            'MaskDic' - Mask dictionary handle.
%                        Default is @MASK.def_bitmask_pipeline.
%            'ExecField' - SIM object field on which to execute the hole
%                        search. Default is 'Im'.
%            'BitType' - Mask bit type. Default is 'uint32'.
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=flag_holes(Sim);
% Reliable: 2

BackField    = SIM.BackField;
ErrField     = SIM.ErrField;

DefV.Nsigma              = 5;
DefV.Bit_Hole            = 'Bit_Hole';
DefV.MaskDic             = @MASK.def_bitmask_pipeline;
DefV.ExecField           = SIM.ImageField;
DefV.BitType             = 'uint32';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (nargout>1)
    SimFlag = SIM(size(Sim));
end

Nsim = numel(Sim);
for Isim=1:1:Nsim
    Tmp = (Sim(Isim).(InPar.ExecField) - Sim(Isim).(BackField))./Sim(Isim).(ErrField) < -InPar.Nsigma;

    Sim(Isim) = bitmask_set(Sim(Isim),Tmp,InPar.Bit_Hole,InPar.BitType);

    if (nargout>1)
        %Sim(Isim),
        SimFlag(Isim).(InPar.ExecField) = Tmp;
    end
end