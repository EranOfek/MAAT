function [Head,UI]=unique(Head)
%--------------------------------------------------------------------------
% unique function                                            class/@Header
% Description: Given an Header object return a new Header
%              with unique keywords names.
% Input  : - Header class object.
% Output : - Header object with unique keywords.
%          - Indices of the unique keywords in the input header.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Head,UI]=unique(Head)
% Reliable: 
%--------------------------------------------------------------------------


HeaderField = 'Header';

% remove empty keywords
Flag = ~Util.cell.isempty_cell(Head.(HeaderField)(:,1));
Head.(HeaderField) = Head.(HeaderField)(Flag,:);

[~,UI]     = unique(Head.(HeaderField)(:,1),'stable');
Head.(HeaderField) = Head.(HeaderField)(UI,:);


