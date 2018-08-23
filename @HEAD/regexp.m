function Out=regexp(Head,Col,varargin)
%--------------------------------------------------------------------------
% regexp function                                              class/@HEAD
% Description: Execute regexp.m on one of the columns of an HEAD object
%              with a single element.
% Input  : - An a single element HEAD object.
%          - Column index on which to execute regexp.m (1|2|3).
%          * Additional arguments to pass to regexp.m, ...,EXPRESSION,...
% Output : - The output of regexp.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: regexp(Head,1,'A_\d+_\d+')
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField    = 'Header';

if (numel(Head)>1),
    error('HEAD/regexp works on a single element HEAD object');
end


Out = regexp(Head.(HeaderField)(:,Col),varargin{:});
