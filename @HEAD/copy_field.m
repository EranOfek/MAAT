function Hout=copy_field(H,FromField,ToField,Hout)
%--------------------------------------------------------------------------
% copy_field function                                          class/@HEAD
% Description: Given an HEAD/AstCat/SIM object copy the content of one
%              field to another field.
% Input  : - An HEAD/AstCat/SIM object from which to copy.
%          - Field name (string) from which to copy.
%          - Field name (string) to copy to.
%          - An HEAD/AstCat/SIM object to which to copy the field.
%            If not provided, then will be the same as the first
%            input argument.
% Output : - An HEAD/AstCat/SIM object with the copied fields.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S1=copy_field(S,'Im','BackIm');
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==3),
    Hout = H;
elseif (nargin==4),
    % Hout is provided by user
else
    error('Illegal number of input arguments: (H,FromField,ToField,[Hout])');
end
    
Nh = numel(H);
for Ih=1:1:Nh,
    Hout(Ih).(ToField) = Hout(Ih).(FromField);
end
    