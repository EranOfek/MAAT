function Data=read_formatted(File,FormatCell,varargin)
% Read text with column position format.
% Package: Util.IO
% Description: Read text/data file with constant format (i.e., columns
%              at specfied locations). The start and end
%              column for each entry should be specified.
% Input  : - File name (string).
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
% See also:  read_str_formatted.m
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jan 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=Util.IO.read_formatted('file.txt',{1,5,'%d'; 7,15,'%f'});
% Reliable: 2
%--------------------------------------------------------------------------


DefV.Skip     = 0;
DefV.Comment  = '%s';

%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);



%Nl   = wc(File,'l');
Ncol = size(FormatCell,1);
%Data = cell(Nl,Ncol);
Data = cell(1,Ncol);

FID = fopen(File,'r');

for Iskip=1:1:InPar.Skip
   fgetl(FID);
end
I   = 0;
while (feof(FID)==0)
   Line = fgetl(FID);
   if (isempty(Line))
      % ignore line
   else
      if (any(strcmp(Line(1),InPar.Comment)))
         % ignore line
      else
         I = I + 1;

         for Icol=1:1:Ncol
            switch FormatCell{Icol,3}
                case {'%s','%c'}
                    Data{I,Icol} = sscanf(Line(FormatCell{Icol,1}:FormatCell{Icol,2}),FormatCell{Icol,3});
                otherwise
                    Data{I,Icol} = Util.string.str2num_nan(Line(FormatCell{Icol,1}:FormatCell{Icol,2}));
            end
         end
      end
   end
end
fclose(FID);

Data = Data(1:I,:);

