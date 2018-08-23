function Head=delete_keytemp(Head,Templates)
% Delete header keywords matching a regular expression.
% Package: @HEAD
% Description: Delete header keywords matching a regular expression.
% Input  : - An HEAD object.
%          - A string or a cell array of strings containing regular
%            expressions. Header keywords corresponding to the the
%            regular expression will be deleted from the header.
% Output : - An HEAD object with the deleted keyword lines.
% See also: delete_key.m, delete_wcs.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Head=delete_keytemp(Head,{'A_\d+_\d+'});
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField    = 'Header';

if (~iscell(Templates))
    Templates = {Templates};
end

Col = 1;

Nh = numel(Head);
Nt = numel(Templates);
for Ih=1:1:Nh
    
    % delete templates
    Flag = false(size(Head(Ih).(HeaderField),1),1);
    for It=1:1:Nt
        Flag = Flag | ~Util.cell.isempty_cell(regexp(Head(Ih).(HeaderField)(:,Col),Templates{It}));
    end
    Head(Ih).(HeaderField) = Head(Ih).(HeaderField)(~Flag,:);   
end
