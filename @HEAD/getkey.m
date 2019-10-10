function [Val,Comment]=getkey(Head,Keyword,SpaceDel,Conv2Num)
% Get a single keyword value from an image header object.
% Package: @HEAD
% Description: Get a single keyword value from an image header object.
% Input  : - Header class object.
%          - Keyword exact name.
%          - Delete spaces from string values {true|false}.
%            Default is true.
%          - A logical flag indicating if to attempt to convert the value
%            to a number (if possible). Default is false.
% Output : - Keyword value. Return a NaN if not exist.
%          - Keyword comment.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Val,Comment]=getkey(Head,'EXPTIME')
% Reliable: 2
%--------------------------------------------------------------------------

Def.SpaceDel = true;
Def.Conv2Num = false;
if (nargin==2)
    SpaceDel = Def.SpaceDel;
    Conv2Num = Def.Conv2Num;
elseif (nargin==3)
    Conv2Num = Def.Conv2Num;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments: getkey(Head,Keyword,[SpaceDel,Conv2Num])');
end

        
HeaderField = 'Header';
Nh      = numel(Head);
Val     = cell(size(Head));
Comment = cell(size(Head));
for Ih=1:1:Nh
    if (~isempty(Head(Ih).(HeaderField)))
        Flag = strcmp(Head(Ih).(HeaderField)(:,1),Keyword);
        if (any(Flag))
            Val{Ih}     = Head(Ih).(HeaderField){Flag,2};
            Comment{Ih} = Head(Ih).(HeaderField){Flag,3};
        else
            Val{Ih}     = NaN;
            Comment{Ih} = NaN;
        end
    else
        Val{Ih}     = NaN;
        Comment{Ih} = NaN;
    end
    if (ischar(Val{Ih}) && SpaceDel)
        Val{Ih} = Util.string.spacedel(Val{Ih});
    end
    
    % attempt to convert to number
    if (Conv2Num)
        Tmp = str2double(Val{Ih});
        if (~isnan(Tmp))
            Val{Ih} = Tmp;
        end
    end
end
            
