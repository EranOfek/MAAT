function Flag=isastcat_col(AstC,Col)
% Check if a string is a valid column name in catalog.
% Package: @AstCat
% Description: Check if a string is a valid column name in catalog.
% Input  : - AstCat object.
%          - A string containing a single column name.
% Output : - An array of logical flags indicating if the column exist in
%            each one of the AstCat object elements.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Flag=isastcat_col(AstC,'XWIN_IMAGE');
% Reliable: 2
%--------------------------------------------------------------------------

ColCellField = 'ColCell';


Nc   = numel(AstC);
Flag = false(size(AstC));

for Ic=1:1:Nc,
    Flag(Ic) = any(strcmp(AstC(Ic).(ColCellField),Col));
end