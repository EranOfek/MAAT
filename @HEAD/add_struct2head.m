function Head=add_struct2head(Head,St)
% Add a structure array into an HEAD object.
% Package: @HEAD
% Description: Add a structure array into an HEAD object.
%              The keyword and values in each element of the structure
%              array will replace or added the keyword in the corresponding
%              HEAD element.
% Input  : - An HEAD object.
%          - A structure array.
% Output : - An HEAD object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Head=add_struct2head(S(1),A)
% Reliable: 2
%--------------------------------------------------------------------------


Nh   = numel(Head);
Nst  = numel(St);

if ~(Nh==Nst || Nh==1 || Nst==1)
    error('Header object and structure have incompatible sizes');
end
Keys      = fieldnames(St);
Nkeys     = numel(Keys);

N = max(Nst,Nh);

for I=1:1:N
    Ist = min(I,Nst);
    Ih  = min(I,Nh);
    
    Vals = struct2cell(St(Ist));
    
    SubHead = [Keys,Vals,cell(Nkeys,1)];
    Head    = replace_key(Head,SubHead);
end

