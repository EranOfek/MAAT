function Head=regexprep(Head,Col,varargin)
%--------------------------------------------------------------------------
% regexprep function                                           class/@HEAD
% Description: Execute regexprep.m on one of the columns of an HEAD
%              object. Unlike HEAD/regexp.m this function is not limited
%              to single element HEAD object.
% Input  : - An HEAD object.
%          - Column index on which to execute regexp.m (1|2|3).
%          * Additional arguments to pass to regexpexp.m,
%            ...,EXPRESSION,REPLACE,...
% Output : - The output of regexp.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: regexp(Head,1,'A_\d+_\d+','BB')
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField    = 'Header';

Nh = numel(Head);
for Ih=1:1:Nh,
    Head(Ih).(HeaderField)(:,Col) = regexprep(Head(Ih).(HeaderField)(:,Col),varargin{:});
end