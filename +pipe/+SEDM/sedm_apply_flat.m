function SI=sedm_apply_flat(SI,FlatSI,varargin)
%--------------------------------------------------------------------------
% sedm_apply_flat function                                            SEDM
% Description: 
% Input  : - 
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

DefV.SpecField     = 'SpexSpecFit';
DefV.FlatSpecField = 'FlatSpexSpecFit';
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

Nseg = numel(SI);
for Iseg=1:1:Nseg,
    SI(Iseg).(InPar.FlatSpecField) = SI(Iseg).(InPar.SpecField)./FlatSI.NormFlat(Iseg);
    SI(Iseg).NormFlat = FlatSI.NormFlat(Iseg);
end

   