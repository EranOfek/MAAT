function OutC=merge(AstC)
%--------------------------------------------------------------------------
% merge function                                             class/@AstCat
% Description: Merge multiple AstCat elements which have the same column
%              structure.
% Input  : - AstCat object with multiple elements and each element have the
%            same number of columns.
% Output : - An AstCat object with the merged rows of the input AstCat
%            object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A=AstCat; A(1).Cat=ones(5,2); A(2).Cat=ones(2,2).*2; A(3).Cat=ones(3,2).*3; 
%          B = merge(A)
% Reliable: 2
%--------------------------------------------------------------------------

Rem0col  = true;

Ncat = numel(AstC);
OutC = AstCat;
OutC.ColCell  = AstC(1).ColCell;
OutC.Col      = AstC(1).Col;
OutC.ColUnits = AstC(1).ColUnits;

[SizeRow,SizeCol] = sizecat(AstC);
Nrow = sum(SizeRow);
Ncol = SizeCol(1);
if (length(unique(SizeCol(SizeCol~=0)  ))~=1)
    error('AstCat elements have different number of columns');
end


if (Rem0col)
    Flag    = SizeCol~=0 & SizeRow~=0;
    AstC    = AstC(Flag);
    SizeCol = SizeCol(Flag);
    SizeRow = SizeRow(Flag);
end



Line     = 0;
Ncat = numel(AstC);
FirstFound = false;
for Icat=1:1:Ncat
   
    if ~isempty(AstC(Icat).Cat)
        % initialization
        if ~FirstFound
            if istable(AstC(Icat).Cat)
                EmptyCell = cell(1,Ncol);
                OutC.Cat = table(EmptyCell{:});
                OutC.Cat.Properties =  AstC(Icat).Cat.Properties;
            else
                OutC.Cat = zeros(Nrow,Ncol);
            end
            FirstFound = true;
        end
        
        
        
        if istable(AstC(Icat).Cat)
            OutC.Cat = [OutC.Cat; AstC(Icat).Cat];
        else
            OutC.Cat((1:SizeRow(Icat))+Line,:) = AstC(Icat).Cat;
        end
        Line = SizeRow(Icat)+Line;
    end
end

