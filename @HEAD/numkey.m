function Nkey=numkey(Head)
% Return the number of keywords (lines) in an header.
% Package: @HEAD
% Description: Return the number of keywords (lines) in an header.
% Input  : - An Header class object.
% Output : - The number of lines (keywords) in the header.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


Nh = numel(Head);
Nkey = zeros(size(Head));
for Ih=1:1:Nh
    Nkey(Ih) = size(Head.Header,1);
end
