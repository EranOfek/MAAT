function Res=wc(FileName,Option)
% Apply the Unix wc (word count) command for a file name.
% Package: Util.files
% Description: Apply the Unix wc (word count) command for a file name in
%              the current directory.
% Input  : - Character array containing file name.
%          - Options:
%            'c' - Return number of bytes only (fast).
%            'm' - Return number of chars only.
%            'l' - Return number of lines only (fast) - default.
%            'L' - Return max line length only.
%            'w' - Return number of words only.
% Output : - Requested number.
% Tested : Matlab 5.3 (on Linux).
%     By : Eran O. Ofek                    Jul 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=wc('wc.m','c');
%--------------------------------------------------------------------------
if (nargin==1),
   Option = 'l';
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input argument');
end

RunStr = sprintf('wc -%c %s',Option,FileName);
[Status,ResStr]=system(RunStr);
Res = sscanf(ResStr,'%f',1);

