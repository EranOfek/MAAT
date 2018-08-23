function [IsVal,KeyVal]=istype(Head,Val,Dictionary)
% Check if TYPE/IMGTYP/OBSTYPE keyword equal some specific value.
% Package: @HEAD
% Description: Check if TYPE/IMGTYP/OBSTYPE keyword equal some specific
%              value.
% Input  : - An HEAD object (or e.g., a SIM image).
%          - Optional image type value, or a cell array of values to check.
%            Values are always strings.
%          - Dictionary of keyword names to retrieve. The first populated
%            keyword will be used.
%            The program check for the value of each one of these
%            keyword values.
%            Default is {'TYPE','IMTYPE','OBSTYPE','IMGTYP','IMGTYPE','IMAGETYP','OBJECT'}.
%            If empty use default.
% Output : - A vector of flags indicating if each one of the HEAD
%            element type is equal to the requested type.
%          - A cell array with the HEAD elements types.
%            Line per HEAD element, column per dictionary keyword.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [IsVal,KeyVal]=istype(Sim,'object')
%          [IsVal,KeyVal]=istype(Sim,'bias')
%          [~,KeyVal]=istype(Sim)
% Reliable: 2
%--------------------------------------------------------------------------

Def.Val        = [];
Def.SubStr     = false;
Def.Dictionary = {'TYPE','IMTYPE','OBSTYPE','IMGTYP','IMGTYPE','IMAGETYP','OBJECT'};
if (nargin==1)
    Val        = Def.Val;
    SubStr     = Def.SubStr;
    Dictionary = Def.Dictionary;
elseif (nargin==2)
    SubStr     = Def.SubStr;
    Dictionary = Def.Dictionary;
elseif (nargin==3)
    Dictionary = Def.Dictionary;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end


if (isempty(Dictionary))
    Dictionary = Def.Dictionary;
end

KeyVal = mgetkey(Head,Dictionary);
if (isempty(Val))
    % do not check - just return KeyVal as is.
    IsVal = [];
else
    if (~iscell(Val))
        Val = {Val};
    end
    IsVal = false(numel(Head),1);
    Nval  = numel(Val);
    for Ival=1:1:Nval
        
        IsVal = IsVal | any(strcmp(KeyVal,Util.string.spacedel(Val{Ival})),2);
    end
end

