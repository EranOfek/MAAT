 function C=mergeStruct(A,B)
    % Merge two structures (unique fields)
    % Input  : - First structure
    %          - Second structure
    % Output : - Structure with merged fields.
    % Eran Ofek      Mar 2021
    [AllFields,I] = unique([fieldnames(A);fieldnames(B)]);
    AllVal        = [struct2cell(A);struct2cell(B)];
    AllVal        = AllVal(I);
    C             = cell2struct(AllVal,AllFields);
end        

