function Sim=flag2nan(Sim,varargin)
%--------------------------------------------------------------------------
% flag2nan function                                             class/@SIM
% Description: Given a SIM image and a flag map (true/false logical) or
%              a bit mask, set all the flagged pixels into NaN.
% Input  : - SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - Cell array of fields, or a field name, in the
%                          input SIM array in which to execute this
%                          function. Default is 'Im'.
%            'FlagMap'   - A User supplied Flag Map which is either a
%                          SIM in which the Im field contains logical
%                          flags, a matrix of logical flags, or a matrix
%                          of bit masks. If empty, then will attempt to use
%                          the .Mask field in the input SIM.
%                          Default is empty.
%            'BitMask'   - A bit mask. If any of this bit mask are
%                          satisfied by the bitmask of the input image
%                          than the pixel will be set to NaN.
%                          Default is 0.
%            'BitType'   - Bit type. If the BitMask in the SIM is numeric
%                          then this parameter is used to specify if the
%                          bit mask is a vector of indices ('index') or a
%                          bit mask value ('value'). Default is 'value'.
% Output : - A SIM object in which the flagged pixels are set to NaN.
% See also: bitmask_find.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=flag2nan(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField   = 'Im';


DefV.ExecField          = {ImageField};
DefV.FlagMap            = [];
DefV.BitMask            = 0;
DefV.BitType            = 'value';  % 'value'|'index' - see bitmask_find.m

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

if (isempty(InPar.FlagMap))
    % User didn't supply FlagMap
    % attempt to use Mask image
    SimFlag = bitmask_find(Sim,InPar.BitMask,'bitand',InPar.BitType);
    % SimFlag is a SIM in which the Im field contains logical flag for pixels
    % that satisfy the bitmask
else
    % User supplied FlagMap which is either a:
    % SIM in which the Im field contains logical flags
    % A matrix of logical flags
    % A matrix of bitmask
    if (SIM.issim(InPar.FlagMap))
        SimFlag = InPar.FlagMap;
    elseif (islogical(InPar.FlagMap))
        SimFlag = SIM;
        SimFlag.(ImageField) = InPar.FlagMap;
    else
        % assume FlagMap is a bit mask
        SimFlag = SIM;
        SimFlag.(ImageField) = bitand(InPar.FlagMap,InPar.BitMask)>0;
    end
end
Nsf =numel(SimFlag);
            

Nsim = numel(Sim);
for Isim=1:1:Nsim
    % for each SIM image
    for If=1:1:Nf
        % for each SIM field
        Isf = min(Isim,Nsf);
        Sim(Isim).(InPar.ExecField{If})(SimFlag(Isf).(ImageField)) = NaN;
    end
end
