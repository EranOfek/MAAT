function Frac=hour_str2frac(String)
% Convert hour string to fraction of day
% Package: statclass/@convert
% Description: Convert a string or cell array of strings containing
%              the hour in format HH:MM:SS.frac' to fraction of day.
% Input  : - String or cell array of strings containing the hour in
%            format HH:MM:SS.frac'.
% Output : - Vector of fraction of day for each hour in the input.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Frac]=convert.hour_str2frac('10:10:10.111')
%          [Frac]=convert.hour_str2frac({'10:10:10.111','11:11:11.12'})
% Reliable: 2
%--------------------------------------------------------------------------

if (ischar(String)),
    String = {String};
end

Time = regexp(String,':','split');
FunH = @(C) (str2double(C{1}).*3600 + str2double(C{2}).*60 + str2double(C{3}))./86400;
Frac = cellfun(FunH,Time);

