function Cell=struct2keyvalcell(Struct)
% Structure field name and content to cell array of key,val pairs
% Package: Util.struct
% Description: Given a structure convert the field names and their values
%              to a cell array of ...,key,val,... arguments.
% Input  : - Structure.
% Output : - A cell array.
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Sep 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S.In=1; S.Out='a'; Cell=Util.struct.struct2keyvalcell(S);
% Reliable: 2
%--------------------------------------------------------------------------

if (isempty(Struct))
    Cell = {};
else
    Key  = fieldnames(Struct);
    Val  = struct2cell(Struct);
    N    = length(Key);
    Cell = cell(1,2.*N);
    for I=1:1:N
        Cell{I.*2-1} = Key{I};
        Cell{I.*2}   = Val{I};
    end
end