function Index=find_strcmpi(varargin)
%--------------------------------------------------------------------------
% find_strcmpi function                                            General
% Description: find(strcmpi(varargin{:})) function. I.e., like strcmpi.m
%              but returning the indices of true.
% Input  : * See strcmpi.m for options.
% Output : - Indices.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: find_strcmpi({'a','b'},'a');
% Reliable: 1
%--------------------------------------------------------------------------


Index = find(strcmpi(varargin{:}));
