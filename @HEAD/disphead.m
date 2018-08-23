function HeadStr=disphead(Head)
% Print headers in HEAD object to screen.
% Package: @HEAD
% Description: Given an HEAD object (or e.g., a SIM image), print all
%              headers to screen.
% Input  : - An HEAD object.
% Output : - A string containing all the headers ready to display.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: disphead(Sim)
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField = 'Header';

Nh = numel(Head);

HeadStr = '';

for Ih=1:1:Nh
    % for each image/header
    HeaderCell = Head(Ih).(HeaderField);
    
    HeadStr = sprintf('%s \n\n Header number %d\n',HeadStr,Ih);
    N = size(HeaderCell,1);
    for I=1:1:N
       if (isempty(HeaderCell{I,1}) && isempty(HeaderCell{I,2}) && isempty(HeaderCell{I,3})),
          % empty line - supress
       else
          if (ischar(HeaderCell{I,2}))
              HeadStr = sprintf('%s%8s%20s%s\n',HeadStr,HeaderCell{I,1},HeaderCell{I,2},HeaderCell{I,3});
          else
              HeadStr = sprintf('%s%8s%20.10f%s\n',HeadStr,HeaderCell{I,1},HeaderCell{I,2},HeaderCell{I,3});
          end
       end
    end
end
