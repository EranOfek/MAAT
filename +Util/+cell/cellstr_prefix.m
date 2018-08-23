function Array=cellstr_prefix(Array,Prefix)
% Add a prefix to all the strings in a char array.
% Package: Util.cell
% Description: Add a prefix to all the strings in a char array.
% Input  : - A char array.
%          - A string prefix.
% Output : - A new char array with prefix.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Array=Util.cell.cellstr_prefix({'a','b','c'},'A_')
% Reliable: 2
%--------------------------------------------------------------------------

if (~iscell(Array)),
    Array = {Array};
end
N = numel(Array);
for I=1:1:N,
    Array{I} = sprintf('%s%s',Prefix,Array{I});
end
