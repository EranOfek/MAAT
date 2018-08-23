function [Key,Val,Comment]=searchkey(Head,String)
%--------------------------------------------------------------------------
% searchkey function                                           class/@HEAD
% Description: Given an Header object, search for keywords name that
%              contains a specific substring.
% Input  : - Header object.
%          - String to search.
% Output : - Candidates header keyword names.
% See also: getkey.m, mgetkey.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Key,Val,Comment]=searchkey(Head,'EXP')
% Reliable: 2
%--------------------------------------------------------------------------

if (numel(Head)>1),
    error('searchkey works on a single element HEAD object');
end

HeaderField = 'Header';
if (~isempty(Head.(HeaderField))),
    F = regexp(Head.(HeaderField)(:,1),String,'match');
    Key = Head.(HeaderField)(~Util.cell.isempty_cell(F),1);
    disp(Key')
    [Val,Comment] = mgetkey(Head,Key);
else
    Key     = {};
    Val     = {};
    Comment = {};

end
