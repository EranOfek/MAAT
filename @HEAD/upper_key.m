function Head=upper_key(Head)
% Convert the key column (first column) of an HEAD object to upper case.
% Package: @HEAD
% Description: Convert the key column (first column) of an HEAD object
%              to upper case.
% Input  : - An Head object.
% Output : - An Head object in which the keyword are in lower case.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Head=upper_key(Head)
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField = 'Header';

Nh = numel(Head);
for Ih=1:1:Nh
    IndChar = cellfun(@ischar,Head(Ih).Header(:,1));
    Head(Ih).(HeaderField)(IndChar,1) = upper(Head(Ih).(HeaderField)(IndChar,1));
end

