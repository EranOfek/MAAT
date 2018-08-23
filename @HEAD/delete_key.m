function Head=delete_key(Head,varargin)
% Delete lines with specific keyword names from an Header object.
% Package: @HEAD
% Description: Delete lines with specific keyword names from
%              an Header object.
% Input  : - An Header object.
%          * Arbitary number of keyword names to delete,
%            or a cell array of keywords to delete.
% Output : - AN Header object with the deleted lines.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Head=delete_key(Head,'EXPTIME');
%          Head=delete_key(Head,{'EXPTIME','AEXPTIME'});
%          Head=delete_key(Head,'AA','BB');
% Reliable: 
%--------------------------------------------------------------------------


HeaderField = 'Header';

if (numel(varargin)==1)
    if (iscell(varargin{1}))
        Key = varargin{1};
    else
        Key = {varargin};
    end
else
    Key = varargin;
end
Nkey = numel(Key);        

Nh = numel(Head);
for Ih=1:1:Nh
    for Ikey=1:1:Nkey
        Flag = strcmpi(Head(Ih).(HeaderField)(:,1),Key{Ikey});
        Head(Ih).(HeaderField) = Head(Ih).(HeaderField)(~Flag,:);
    end
end



