function [Val,Comment,SubHead]=mgetkey(Head,Keyword,SpaceDel,Conv2Num)
% Get imgae header keyword values for multiple keywords.
% Package: @HEAD
% Description: Given an Header object, get imgae header keyword values for
%              multiple keywords. Also get a sub header with the requested
%              key/val/comment.
% Input  : - Header class object, or e.g., SIM object.
%          - Cell array of keyword exact name.
%          - Delete spaces from string values {true|false}.
%            Default is true.
%          - A logical flag indicating if to attempt to convert the value
%            to a number (if possible). Default is false.
% Output : - Cell array of keywords value. Line per HEAD element,
%            column per keyword. Return NaN for non existing keywords.
%          - Cell array of keywords comment.
%          - A sub header containing only the requested key/val/comment.
% See also: getkey.m, searchkey.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Val,Comment,SubHead]=mgetkey(Head,{'EXPTIME','AEXPTIME'});
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
    error('Illegal number of input arguments: mgetkey(Head,Keyword,[SpaceDel,Conv2Num])');
end

HeaderField = 'Header';


if (~iscell(Keyword))
    Keyword = {Keyword};
end

Nh       = numel(Head);
Nkey     = numel(Keyword);
Val      = cell(Nh,Nkey);
Comment  = cell(Nh,Nkey);

SubHaed = HEAD;
for Ih=1:1:Nh
    Isub = 0;
    if (~isempty(Head(Ih).(HeaderField)))
        %Val  = cell(Nkey,1);
        %Comment = cell(Nkey,1);
        for Ikey=1:1:Nkey
            Flag = strcmp(Head(Ih).(HeaderField)(:,1),Keyword{Ikey});
            if (any(Flag))
                Val{Ih,Ikey}     = Head(Ih).(HeaderField){Flag,2};
                Comment{Ih,Ikey} = Head(Ih).(HeaderField){Flag,3};
                if (nargout>2)
                    Isub = Isub + 1;
                    SubHead(Ih).(HeaderField)(Isub,:) = Head(Ih).(HeaderField)(Flag,:);
                end
            else
                Val{Ih,Ikey}     = NaN;
                Comment{Ih,Ikey} = NaN;
            end
            if (ischar(Val{Ih,Ikey}) && SpaceDel)
                Val{Ih,Ikey} = Util.string.spacedel(Val{Ih,Ikey});
            end

            % attempt to convert to number
            if (Conv2Num)
                Tmp = str2double(Val{Ih,Ikey});
                if (~isnan(Tmp))
                    Val{Ih,Ikey} = Tmp;
                end
            end

        end
    else
        Val(Ih,:)  = num2cell(nan(Nkey,1));
        Comment(Ih,:) = num2cell(nan(Nkey,1));
    end
    
end

