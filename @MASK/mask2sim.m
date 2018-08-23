function Sim=mask2sim(Mask,Sim,MaskDic)
% Copy a MASK object into a new or existing SIM object.
% Package: @SIM
% Description: Copy a MASK object into a new or existing SIM object.
% Input  : - A MASK object.
%          - A SIM object. If empty, then will copy the MASK into a new
%            SIM object. Default is empty.
%          - MaskDic to use if the .MaskDic field is not populated.
%            Default is @MASK.def_bitmask_pipeline.
% Output : - A SIM object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=mask2sim(Mask,Sim);
% Reliable: 2
%--------------------------------------------------------------------------

MaskField     = 'Mask';
MaskDicField  = 'MaskDic';

Def.Sim     = [];
Def.MaskDic = @MASK.def_bitmask_pipeline;
if (nargin==1)
    Sim     = Def.Sim;
    MaskDic = Def.MaskDic;
elseif (nargin==2)
    MaskDic = Def.MaskDic;
elseif (nargin==3)
    % do nothing
else
    error('Illegal number of input arguments: [Sim]=mask2sim(Mask,[Sim,MaskDic])');
end


Nmask = numel(Mask);
if (isempty(Sim))
    Sim = SIM(size(Mask));
end
Nsim = numel(Sim);
if (Nmask~=Nsim)
    error('SIM and MASK size must be the same');
end
for Imask=1:1:Nmask
    Sim(Imask).(MaskField)    = Mask(Imask).(MaskField);
    if (isempty(Mask(Imask).(MaskDicField)))
        Sim(Imask).(MaskDicField) = MaskDic;
    else
        Sim(Imask).(MaskDicField) = Mask(Imask).(MaskDicField);
    end
end
    