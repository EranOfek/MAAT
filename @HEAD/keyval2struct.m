function [StVal,StCom]=keyval2struct(Head,Keys,OutKeys,SpaceDel)
%--------------------------------------------------------------------------
% keyval2struct function                                       class/@HEAD
% Description: Get the value of header keywords from an HEAD object
%              and return the result in a structure format.
% Input  : - An HEAD object (or e.g., a SIM image).
%          - Cell array of header keywords to extract from headers.
%          - The field names in which to store the output keyword values.
%            Default is to use the names of the requested keys.
%            This is useful if the keyword name contains an illegal
%            character for a matlab field name (e.g., '-').
%          - Delete spaces from string values {true|false}.
%            Default is true.
% Output : - Structure array (element per HEAD element) of keyword values.
%          - Structure array (element per HEAD element) of comment values.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [StVal,StCom]=keyval2struct(S,{'EXPTIME','FILTER'});
%          [StVal,StCom]=keyval2struct(S,{'EXPTIME','FILTER','SEEING','DATE-OBS'},{'A','B','C','D'})
% Reliable: 2
%--------------------------------------------------------------------------

Def.OutKeys  = Keys;
Def.SpaceDel = true;
if (nargin==2),
    OutKeys  = Def.OutKeys;
    SpaceDel = Def.SpaceDel;
elseif (nargin==3),
    SpaceDel = Def.SpaceDel;
elseif (nargin==4),
    % do nothing
else
    error('Illegal number of input arguments');
end

if (~iscell(Keys)),
    Keys = {Keys};
end
if (~iscell(OutKeys)),
    OutKeys = {OutKeys};
end

if (isempty(OutKeys)),
    OutKeys = Keys;
end

[Val,Com] = mgetkey(Head,Keys,SpaceDel);

StVal = cell2struct(Val,OutKeys,2);
if (nargout>1),
    StCom = cell2struct(Com,OutKeys,2);
end