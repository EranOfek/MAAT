function Num=str2double_check(String)
%--------------------------------------------------------------------------
% str2double_check function                                        General
% Description: Convert strings to doubles, but unlike str2double if the
%              input is a number will return the number (instead of NaN).
% Input  : - String or cell array of strings.
% Output : - Strings converted to numbers.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Num=str2double_check('33');
%          Num=str2double_check({'33',44});
% Reliable: 2
%--------------------------------------------------------------------------

if (iscell(String)),
    
    I = find(cellfun(@isnumeric,String));
    Num = str2double(String);
    if (~isempty(I)),
        Num(I) = cell2mat(String(I));
    end
else
    if (ischar(String)),
        Num = str2double(String);
    elseif (isnumeric(String)),
        Num = String;
    else
        error('Unknwon input option');
    end
end
