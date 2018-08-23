function [IsKeyVal,IsEqual]=iskeyval(Head,varargin)
%--------------------------------------------------------------------------
% iskeyval function                                            class/@HEAD
% Description: Check if header keywords in an HEAD object are equal to
%              some values. For example, you can use this to select all
%              images which are 60 s R-band image (see example).
% Input  : - An HEAD object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword is the header keyword name, and value is its
%            requested value. For string do case senstive comparison.
% Output : - A vector of flags indicating, per each HEAD element, if its
%            all its keywords are equal to the requested values.
%          - A matrix of flags (line per HEAD element, column per keyword)
%            indicating if the value of the keyword is equal to the
%            requested value. 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: IsKeyVal=iskeyval(Head,'SLITGRAB','deployed','DICHTRAN','deployed')
%          IsKeyVal=iskeyval(Head,'EXPTIME',60,'FILTER','R')
% Reliable: 2
%--------------------------------------------------------------------------


Narg = numel(varargin);
Keys = varargin(1:2:Narg-1);
Vals = varargin(2:2:Narg);
Nkey = numel(Keys);

KeyVal = mgetkey(Head,Keys);

Nh = numel(Head);
IsEqual = false(Nh,Nkey);
for Ih=1:1:Nh,
    % for each HEAD element
    for Ikey=1:1:Nkey,
        % for each keyword
        HeadKeyVal   = KeyVal{Ih,Ikey};
        RequestedVal = Vals{Ikey};
        if (ischar(HeadKeyVal) && ischar(RequestedVal)),
            % value is string
            IsEqual(Ih,Ikey) = strcmp(Util.string.spacedel(HeadKeyVal),RequestedVal);
        elseif (isnumeric(HeadKeyVal) && isnumeric(RequestedVal)),
            % value is numeric
            IsEqual(Ih,Ikey) = HeadKeyVal==RequestedVal;
        else
            % do nothing - already false
        end
    end
end
% check if condition is satisfied
IsKeyVal = all(IsEqual,2);
