function [Out,FileCell]=for_each_file(FileName,Function,SaveType,varargin)
% Execute a function on a list of files.
% Package: Util.files
% Description: Given a file name containing list of files, load each file
%              into a matrix and execute a function with the loaded
%              matrix as a paramter.
% Input  : - File name containing a list of files to load.
%            Ignore lines begining with # or %.
%          - Function to execute for each file that is being loaded.
%          - SaveType:
%            'c'   - Save the output from the executed function in a cell
%                    array, default.
%            'm'   - Save the output from the executed function in a matrix,
%                    each line per each execuation (assuming the output is
%                    a raw vector.
%          * Additional aribitrary number of argument to pass to the
%            function to execute.
% Output : - Cell array or matrix containing the output of executing
%            the function on each file.
%          - Cell array of all the individual file names in the input file.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Dec 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%-----------------------------------------------------------------------------
if (nargin==2)
   SaveType = 'c';
end

FID = fopen(FileName,'r');

I = 0;
while (feof(FID)==0)
   Line = fgetl(FID);

   if (isempty(strfind(Line(1),'%'))==1 && isempty(strfind(Line(1),'#'))==1)
      I = I + 1;
      FileCell{I} = Line;

      OutF = feval(Function,load(Line),varargin{:});

      switch SaveType
       case 'c'
	  Out{I} = OutF;
       case 'm'
         if (isempty(OutF)==1)
	    OutF = NaN;    % a partial solution!
         end
	 Out(I,:) = OutF;
       otherwise
          error('Unknown SaveType Option');
      end

   else
      % comment - ignore
   end
end
