function Val=bitand_array(Array,Dim)
% Perform a bitand operation along all elements in an array.
% Package: Util.array
% Description: Perform a bitand operation along all elements in an array
%              along a specific dimension.
% Input  : - An array of integers.
%          - Dimension along to perform the bitand operation. Default is 1.
% Output : - The result of the bitand operation.
% See also: sum_bitor.m (the same)
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Val=bitand_array(Array);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
    Dim = 1;
end

C = class(Array);
switch lower(C)
    case {'uint8','int8'}
        Nbit = 8;
    case {'uint16','int16'}
        Nbit = 16;
    case {'uint32','int32'}
        Nbit = 32;
    case {'uint64','int64'}
        Nbit = 64;
    otherwise
        error('Unknown class - only integers are allowed');
end

Val = 0;
for Ibit=1:1:Nbit
    Val = Val + (2.^(Ibit-1)).*all(bitget(Array,Ibit),Dim);
end
    
    
    
    