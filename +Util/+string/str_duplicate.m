function NewStr=str_duplicate(Str,Num,StrLast)
% Duplicate a string multiple times.
% package: Util.string
% Description: Duplicate a string multiple times.
% Input  : - The string to duplicate.
%          - Number of times to duplicate the string.
%          - Additional string to add at the end of the output string.
%            Default is ''.
% Output : - Duplicated string.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Nov 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: NewStr=str_duplicate('%f,',5,'%f\n');
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2),
    StrLast = '';
end

NewStr = '';
for I=1:1:Num,
    NewStr = sprintf('%s%s',NewStr,Str);
end

NewStr = sprintf('%s%s',NewStr,StrLast);
