function Cat=replace_col_name(Cat,varargin)
% Replace column names in AstCat object
% Package: @AstCat
% Description: Replace column names in AstCat object
% Input  : - An AstCat object
%          * Arbitrary number of pairs of arguments:
%            ...,old_col_name,new_col_name,...
% Output : - An AstCat object
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S
% Reliable: 
%--------------------------------------------------------------------------

ColCellField = AstCat.ColCellField;

Nvar  = numel(varargin);
if (Nvar.*0.5)~=floor(Nvar.*0.5)
    error('Must provide pairs of old_key, new_key');
end

Ncat = numel(Cat);
for Icat=1:1:Ncat
    
    for Ivar=1:2:Nvar-1
        Icol = find(strcmp(Cat(Icat).(ColCellField),varargin{Ivar}));
        if (~isempty(Icol))
            Cat(Icat).(ColCellField){Icol} = varargin{Ivar+1};
        end
    end
end

Cat = colcell2col(Cat);

