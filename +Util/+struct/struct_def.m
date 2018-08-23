function Struct=struct_def(Fields,N,M)
% Define a structure array of a specific size with fields.
% Package: Util.struct
% Description: Define a structure array of a specific size with fields
%              specified in a cell array of names.
% Input  : - A cell array of field names.
%          - Number of lines in array, or [lines, rows]. Default is 1.
%          - Number of Rows in array. Default is 1.
% Output : - A structure array.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Struct=Util.struct.struct_def({'Field1','Field2'},2,1);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
    N = 1;
    M = 1;  
elseif (nargin==2)
    if (numel(N)>1)
        M = N(2);
        N = N(1);
    else
        M = 1;
    end      
else
    % do nothing
end

Nf   = numel(Fields);
Pars = cell(1,Nf.*2);
Ind = 0;
for If=1:2:2.*Nf
    Ind = Ind + 1;
    Pars{If} = Fields{Ind};
    Pars{If+1} = cell(N,M);
end

Struct = struct(Pars{:});

