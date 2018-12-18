function TrimStr=spacedel(Str,Numeric)
% recursively delete all spaces from a string.
% Package: Util.string
% Description: Given a string, recursively delete all spaces.
% Input  : - A string, or a cell array of strings.
% Output : - A trimmed string.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Apr 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: spacetrim.m
% Example: spacedel('a   a ');
%          spacedel({'a   a ';' bbbaj   a'});
% Reliable: 2
%---------------------------------------------------------------------------

TrimStr = regexprep(Str,' ','');
if (iscell(TrimStr))
   If = find(Util.cell.isempty_cell(strfind(TrimStr,' '))==0);
else
   If = strfind(TrimStr,' ');
end
if (~isempty(If))
   TrimStr=spacetrim(TrimStr);
end


