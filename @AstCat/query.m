function [AstC,Flag]=query(AstC,QueryString,varargin)
% Query a catalog using operators on its columns
% Package: @AstCat
% Description: Query a single element AstCat object by evaluating a string
%            of operators on the catalog column names.
% Input  : - An AstCat object.
%          - A query string. E.g., 'RA>0 & RA<1 & Dec>1' or a numeric
%            operation (e.g., 'RA+Dec').
% Output : - An AstCat object with the selected rows.
%          - A vector of logical flags indicating the selected rows, or a
%            vector of the numerical results.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A = query(TD1,'RA>0 & RA<0.1 & V<10');
% Reliable: 2
%--------------------------------------------------------------------------


ColCellField = AstCat.ColCellField;
ColField     = AstCat.ColField;
CatField     = AstCat.CatField;

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (numel(AstC)>1)
    error('AstCat object must contain a single element');
end


% for table:  'RA>1' -> 'AstC.(CatField).(AstC.Col.RA)'
% for matrix: 'RA>1' -> 'AstC.(CatField)(:,AstC.Col.RA)'

ColStrLength = cellfun(@numel,AstC.(ColCellField));
Ncol   = numel(AstC.(ColCellField));
% need to go over column names from longest to shortest
[~,SI] = sort(ColStrLength,'descend');

QueryString = sprintf(' %s',QueryString);
if istable(AstC.(CatField))
    % catalog is stored as table
    for Icol=1:1:Ncol
        ColName = AstC.(ColCellField){SI(Icol)};
        %QueryString = strrep(QueryString,ColName,sprintf('AstC.(CatField).%s',ColName));
        %QueryString = regexprep(QueryString,sprintf('[^.]%s',ColName),sprintf(' AstC.(CatField).%s',ColName));
        QueryString = regexprep(QueryString,sprintf('(?<=[^.])%s',ColName),sprintf(' AstC.(CatField).%s',ColName));
    end
else
    % catalog is stored as array
    for Icol=1:1:Ncol
        ColName = AstC.(ColCellField){SI(Icol)};
        %QueryString = strrep(QueryString,ColName,sprintf('AstC.(CatField)(:,AstC.(ColField).%s',ColName));
        %QueryString = regexprep(QueryString,sprintf('[^.]%s',ColName),sprintf(' AstC.(CatField)(:,AstC.(ColField).%s)',ColName));
        QueryString = regexprep(QueryString,sprintf('(?<=[^.])%s',ColName),sprintf(' AstC.(CatField)(:,AstC.(ColField).%s)',ColName));
    end
end
Flag = eval(QueryString);

if (islogical(Flag))
    AstC.(CatField) = AstC.(CatField)(Flag,:);
end
    


