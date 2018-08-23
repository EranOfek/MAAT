function Head=update_key(Head,varargin)
%--------------------------------------------------------------------------
% update_key function                                          class/@HEAD
% Description: Update the value and comment of a keyword in an Header
%              object. If keyword doesn't exist then do nothing.
% Input  : - An Header object.
%          * Either an arbitrary number of triplets:
%            ...,key,val,comment,...;
%            or an arbitrary number of cell arrays, each cell array
%            containing 2 or 3 columns of {key,val,comment}, where the
%            comment field is optional;
%            or a single element HEAD object.
% Output : - An Header object in which the content of the keyword was
%            updated.
% See also: replace_key.m, delete_key.m, add_key.m, getkey.m, mgetkey.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Head=replace_key(Head,'EXPTIME',60,'Exp. time [s]','NAXIS',3,'');
%          Head=replace_key(Head,{'EXPTIME',60;'NAXIS',2});
%          H = replace_key(S,S(1));
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField = 'Header';

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
            if (all(~Flag))
                %warning('Keyword %s was not found in header',Key);
                % do nothing
                
            else
                Head(Ih).(HeaderField){Flag,1} = Key;
                Head(Ih).(HeaderField){Flag,2} = Val;
                if (~isempty(Com))
                    Head(Ih).(HeaderField){Flag,3} = Com;
                end
            end
        end
    end
end

    
    
