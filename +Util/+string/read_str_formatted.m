function [Data]=read_str_formatted(Str,FormatCell,varargin)
% Read a string which is formatted by specidc column positions 
% Package: Util.string

% Description: Read text/data string with constant format (i.e., columns
%              at specfied locations). The start and end
%              column for each entry should be specified.
% Input  : - Either a string containing a file content, or a cell array
%            in which each element contains a line from the file.
%          - Cell array containing three columns, specifing the file
%            structure:
%            {start column, end column, format string}.
%            Example: {1,5,'%d'; 7,15,'%f'}.
%          * Arbitrary number of pairs of arguments {...,key,val,...},
%            Where the following keyword are allowed:
%            'skip'   - number of lines to skip, default is 0.
%            'comment'- comment character, or a cell array of comments
%                       characters. Default is '%'.
% Output : - Cell array containing data.
% See also: read_formatted.m
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=Util.string.read_str_formatted(String,{1,5,'%d'; 7,15,'%f'});
% Reliable: 2
%--------------------------------------------------------------------------


DefV.Skip     = 0;
DefV.Comment  = '%s';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (iscell(Str))
    Lines = Str;
    clear Str;
else
    Lines = regexp(Str,'\n','split');
end

Nline = numel(Lines);

Ncol  = size(FormatCell,1);
Data  = cell(Nline,Ncol);

I = 0;
for Il=1:1:Nline
   if (isempty(Lines{Il}) || Il<=InPar.Skip)
      % ignore line
   else
      if (any(strcmp(Lines{Il}(1),InPar.Comment)))
         % ignore line
      else
         I = I + 1;
         for Icol=1:1:Ncol
             if (Icol==Ncol)
                 EndRow = min(FormatCell{Icol,2},length(Lines{Il}));
             else
                 EndRow = FormatCell{Icol,2};
             end
             
            switch FormatCell{Icol,3}
                case {'%s','%c'}
                    Data{I,Icol} = sscanf(Lines{Il}(FormatCell{Icol,1}:EndRow),FormatCell{Icol,3});
                otherwise
                    Data{I,Icol} = Util.string.str2num_nan(Lines{Il}(FormatCell{Icol,1}:EndRow));
            end
         end
      end
   end
end

Data = Data(1:I,:);
