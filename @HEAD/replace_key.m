function Head=replace_key(Head,varargin)
% Replace the value and comment of a keyword in an Header object.
% Package: @HEAD
% Description: Replace the value and comment of a keyword in an Header
%              object. If keyword doesn't exist then it will be added
%              to header.
% Input  : - An Header object.
%          * Either an arbitrary number of triplets:
%            ...,key,val,comment,...;
%            or an arbitrary number of cell arrays, each cell array
%            containing 2 or 3 columns of {key,val,comment}, where the
%            comment field is optional;
%            or a single element HEAD object.
% Output : - An Header object in which the content of the keyword was
%            replaced.
% See also: update_key.m, delete_key.m, add_key.m, getkey.m, mgetkey.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Head=replace_key(Head,'EXPTIME',60,'Exp. time [s]','NAXIS',3,'');
%          Head=replace_key(Head,{'EXPTIME',60;'NAXIS',2});
%          H = replace_key(S,S(1));
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField = HEAD.HeaderField;

if (~iscell(varargin{1}) && ~HEAD.ishead(varargin{1}) )
    % assume varargin contains triplets of key,val,comments
    Ntriplets   = numel(varargin)./3;
    Tmp = reshape(varargin,3,Ntriplets).';
    clear varargin;
    varargin{1} = Tmp;
end

Narg = numel(varargin);


Nh = numel(Head);
for Ih=1:1:Nh
    if (isempty(Head(Ih).(HeaderField)))
        Head(Ih).(HeaderField) = cell(0,3);
    end
    
    for Iarg=1:1:Narg
        if (HEAD.ishead(varargin{Iarg}))
            varargin{Iarg} = varargin{Iarg}.(HeaderField);
        end
        
        
        [Nline,Nrow] = size(varargin{Iarg});
        for Iline=1:1:Nline
            Key = varargin{Iarg}{Iline,1};
            Val = varargin{Iarg}{Iline,2};
            if (Nrow>2)
                Com = varargin{Iarg}{Iline,3};
            else
                Com = '';
            end
            
            % search for key in Header
            Flag = strcmp(Head(Ih).(HeaderField)(:,1),Key);
            %%if (all(~Flag))
            if (all(~Flag)) || strcmpi(Key, 'COMMENT') % Na'ama, 20180907    
                %warning('Keyword %s was not found in header',Key);
                Nline = size(Head(Ih).(HeaderField),1);
                Head(Ih).(HeaderField){Nline+1,1} = Key;
                Head(Ih).(HeaderField){Nline+1,2} = Val;
                Head(Ih).(HeaderField){Nline+1,3} = Com;
                
            else
                if (sum(Flag)>1)
                    error('Two keywords with the same name');
                end
                Head(Ih).(HeaderField){Flag,1} = Key;
                Head(Ih).(HeaderField){Flag,2} = Val;
                Head(Ih).(HeaderField){Flag,3} = Com;
            end
        end
    end
end

    
    
