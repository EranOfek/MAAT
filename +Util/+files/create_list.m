function [ListFileName,ListCell,Exist]=create_list(List,OutputFile,Verify)
% Create a file and a ell array containing a list of files.
% Package: Util.files
% Description: Create a file and a cell array containing a list of files.
%              The list is created from a cell array, or file name
%              with wildcards.
% Input  : - File listing can be one of the followings:
%            (1) A cell vector containing a file name in each cell.
%            (2) A string containing wild cards (i.e., '*' or '?'),
%                in this case, the program will produce a list of
%                all matched files in the current directory.
%                The program uses superdir, so expressions like
%                'lred00[75-82].fits' are legal.
%            (3) A string containing a file name (the file name not
%                necesserly exist).
%            (4) A string begining with '@' containing a file name
%                that contains a list of file (one per line).
%            (5) A matrix or a cell array of matrices.
%                In this case the output in ListCell will be a cell
%                array containing the matrices.
%          - Optional output file name, default is empty (i.e., []).
%            If empty, then create a temporary file name for output.
%            If NaN, than do not create file.
%          - Verify {'y'|'n'} if each one of the files listed in
%            the list exist in the current directory, default is 'n'.
%            If the file doesnot exist then donot print it to the
%            final output list.
% Output : - File name containing the list of files (if verify set to 'y',
%            then this contains only existing files).
%          - A cell array containing the list of files (all files,
%            existing or not).
%          - A vector containing flags indicating if each file in
%            the original list is exist (1) or is missing (0).
%            If the veify option is set to 'n', then return an empty
%            vector.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [F,C]=create_list({'A','B'});
%          [F,C,E]=create_list({'A','B'},[],'y');
%          [~,C]=create_list('lred00[75-82].fits',NaN);
%          [~,C]=create_list('@file',NaN);
%          [~,C]=create_list('A*.fits',NaN);
% Reliable: 1
%--------------------------------------------------------------------------
DefOutputFile = [];
DefVerify      = 'n';

if (nargin==1)
   OutputFile   = DefOutputFile;
   Verify       = DefVerify;
elseif (nargin==2)
   Verify       = DefVerify;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end


if (isempty(List))
   ListFileName = NaN;
   ListCell = {};
   Exist    = [];
else
   IsNumeric = 0;
   IsCell    = 0;
   if (isnumeric(List))
      IsNumeric = 1;
   end
   if (iscell(List))
      if (isnumeric(List{1}))
         IsNumeric = 1;
         IsCell    = 1;
      end
   end
   
   if (IsNumeric)
      % Input is a matrix - return in a cell
      if (IsCell)
         ListCell = List;
      else
         ListCell{1} = List;
      end
      ListFileName = NaN;
      Exist        = 1;
   else
      % input is some sort of a list
      
      if (isempty(OutputFile)==1)
         OutputFile = tempname;
      end
   
      ListFileName = OutputFile;
      
      
      if (ischar(List)==1)
         % search for '@' at the begining
         if (strcmpi(List(1),'@')==1)
            % file containing a list of files...
            List = List(2:end);
            % read files in file name to a cell array
            FID = fopen(List,'r');
            FI = 0;
            while (feof(FID)==0)
                FI = FI + 1;
                ListCell{FI} = fgetl(FID);
            end
            fclose(FID);
         else
            % wild cards or single file
      
            % get parent directory
            DirPos = strfind(List,filesep);
            if (isempty(DirPos)==1)
                DirPath = '';
            else
                DirPath = List(1:DirPos(end));
            end
      
            %LS = dir(List);
            LS = Util.files.superdir(List);
            [ListCell{1:1:length(LS)}] = deal(LS.name);
            if (isempty(LS) && isempty(strfind(List,'*'))==1 && isempty(strfind(List,'?'))==1)
               % single - file not exist
               ListCell{1} = List;
            end
            for Ils=1:1:length(LS)
                ListCell{Ils} = sprintf('%s%s',DirPath,ListCell{Ils});
            end
               
         end
      elseif (iscell(List)==1)
         ListCell = List;
      else
         error('Unknwon format for List');
      end
      
      Nf = length(ListCell);    % number of files in ListCell
      Exist = ones(Nf,1);       % set to all exist
      
      switch lower(Verify)
       case 'y'
          for If=1:1:Nf
             FID_try = fopen(ListCell{If},'r');
             if (FID_try==-1)
                Exist(If) = 0;
             else
                Exist(If) = 1;
                fclose(FID_try);
             end
          end
       case 'n'
          % do nothing
       otherwise
          error('Unknwon Veify option');
      end
      
      % Write output file
      if (~isnan(OutputFile))
         Iex = find(Exist==1);
         SelectedListCell = {ListCell{Iex}};
         Util.IO.fprintf_cell(ListFileName,'%s\n',SelectedListCell);
      end
      
      switch lower(Verify)
       case 'n'
          Exist = [];
       otherwise
          % do nothing
      end
   end
end
