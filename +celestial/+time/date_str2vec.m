function DateVec=date_str2vec(String)
%--------------------------------------------------------------------------
% date_str2vec function                                            General
% Description: Convert a string or a cell array of string containing date
%              and time in the format 'YYYY-MM-DD HH:MM:SS.frac'
%              or 'YYYY-MM-DD', to a matrix of dates with the following
%              columns [Y M D H M S].
%              OBSOLETE: Use convert.str2date instead.
% Input  : - A string or a cell array of string containing date
%            and time in the format 'YYYY-MM-DD HH:MM:SS.frac'
%            or 'YYYY-MM-DD'
% Output : - A matrix of dates with the following columns [Y M D H M S].
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: DateVec=date_str2vec({'2010-10-11T10:10:10.11','2014-12-12 10:23:59.1'});
% Reliable: 2
%--------------------------------------------------------------------------

if (ischar(String))
    String = {String};
end
if (size(String,1)==1)
    String = String.';
end
Time = regexp(String,'T|:|\s|-','split');
DateVec = cell2mat(cellfun(@str2double,Time,'UniformOutput',false));

