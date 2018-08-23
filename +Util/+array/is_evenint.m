function Flag=is_evenint(Int)
% Check for each integer number in array if even.
% Package: Util.array
% Description: Check for each integer number in array if even.
% Input  : - Array of integers.
% Output : - Array of boolean flags indicating if each elemnt is even.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Apr 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Flag=is_evenint([-1 0 1 2 3]);
% Reliable: 2
%--------------------------------------------------------------------------

Flag = Int.*0.5 == floor(Int.*0.5);