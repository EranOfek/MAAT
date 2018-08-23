function CellLines=strlines2cell(String)
%--------------------------------------------------------------------------
% strlines2cell function                                           General
% Description: Given a string with end of line characters, break the
%              string into a cell array in which each cell contains
%              a line.
% Input  : - A string.
% Output : - A cell array of lines.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: CellLines=strlines2cell(String)
% Reliable: 2
%--------------------------------------------------------------------------

CellLines = regexp(String,'\n','split');
