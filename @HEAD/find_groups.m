function Groups=find_groups(Head,Keys,DelSpace)
%--------------------------------------------------------------------------
% find_groups function                                         class/@HEAD
% Description: Find groups of headers which have keywords with the
%              identical values. Given an HEAD object with multiple
%              elements, and a list of keywords find unique groups
%              defined by the same keyword values. For example, images
%              which were taken with the same filter and same exposure
%              time. This function use the Util.cell.cell_find_groups.m function.
% Input  : - An HEAD object.
%          - Cell array of keywords by which to find groups.
%          - Delete spaces from strings prior to comparison {true|false}.
%            Default is true.
% Output : - Structure of groups. Each structure represent a group of rows
%            in the input cell array which have equal values.
%            The structure contains the following fields:
%            .Conntent  - Cell array of values define the group.
%            .ptr       - The indices of the rows that belong to the
%                         group.
% See also: cell_find_groups.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Groups=find_groups(Sim,{'EXPTIME','FILTER'});
%          Groups=find_groups(lower_key(Sim),{'exptime','filter'});
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3)
    DelSpace = true;
end

% construct a cell array of header keywords
% line per image/header, column per keyword:
Cell = mgetkey(Head,Keys);

% find groups:
Groups = Util.cell.cell_find_groups(Cell,DelSpace);

