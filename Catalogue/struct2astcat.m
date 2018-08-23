function AstC=struct2astcat(St)
%--------------------------------------------------------------------------
% struct2astcat function                                         Catalogue
% Description: Convert a structure containing an Astronomical catalog
%              into an AstCat object.
%              OBSOLETE - USE AstCat.struct2astcat instead
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=struct2astcat(FIRST)
% Reliable: 2
%--------------------------------------------------------------------------

if (numel(St)>1)
    error('AstCat should have a single element');
end



if (isastcat(St))
    % object is already in AstCat format
    AstC = St;
else
    AstC = AstCat;
    AstCatFields = fieldnames(AstC);
    
    %Fields = fieldnames(St);
    Nf     = numel(AstCatFields);
    for If=1:1:Nf
        if (isfield(St,AstCatFields{If}))
            AstC.(AstCatFields{If}) = St.(AstCatFields{If});
        end
    end
end