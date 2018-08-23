function [Cat,Msg]=run_galex_sql(Query,OutputFileName,QueryType,ConcatString)
% Run a GALEX command line SQL quary (OBSOLETE - see VO.MAST)
% Package: VO.GALEX
% Description: Run a GALEX command line SQL quary. Running the java
%              command line tool.
% Input  : - File name containing the SQL query, or alternatively,
%            a string containing the full query,  or alternatively,
%            this can be a 3 elements cell array, in which the 
%            first cell contains : a cell array of attributes in the
%                                  SELECT clause.
%            second cell contains: a cell array of tables/views in the
%                                  FROM clause.
%            third cell contains : a string to appear in the WHERE clause. 
%          - Output file name.
%            If empty (i.e., []), then create a tmp file, which will
%            be deleteded at the end of the process, default is [].
%          - Query type:
%            'file'   - First argument is file name containing query (default).
%            'string' - First argument is the full query string.
%            'cell'   - First argument is a cell array of attributes.
%          - Optional string to concatenate to query (e.g., this could be
%            a variable part of the query), default is empty string (i.e., '').
% Output : - Optional catalog, containing the search result.
%          - Return code: 1 - objects found; 0- No objects have been found
%          - First line in resulted file.
% Instellation: Follow the casjobs.jar instellation instructions in:
%            http://galex.stsci.edu/casjobs/casjobscl.aspx
%            install the casjobs.jar and CasJobs.config files in
%            .../matlab/fun/bin/CasJobs/GALEX/
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Needed : The python script sqlcl.py should be located in ./sqlcl/ directory
%          relative to the location of this script.
% Notes  : Generic SDSS queries can be found in the ./sqlcl/sql/ directory.
% Example: Scheme browser at: http://galex.stsci.edu/GR6/?page=dbinfo
%          Static file query:
%          Query='/scr2/eran/bin/sdssQA/CasJOB/sql/get_1004.sql';
%          [Cat,Msg,FirstLine]=VO.GALEX.run_galex_sql(Query,[],'file');
%          % dynamic file query:
%          Query='/scr2/eran/matlab/sdss/bin/sql/get_obj_field.sql';
%          Con='run=3813 and rerun=41 and camcol=6 and field=60';
%          [Cat,Msg,FirstLine]=VO.GALEX.run_galex_sql(Query,[],'file',Con);
%          % string query:
%          QueryStr = 'select top 2 ra,dec from star';
%          [Cat,Msg,FirstLine]=VO.GALEX.run_galex_sql(QueryStr,[],'string');
%          % Cell query:
%          Q{1} = {'ra','dec','nuv_mag','nuv_magerr'};
%          Q{2} = {'photoobjall'};
%          Q{3} = '(ra between 0 and 0.1) and (dec between 0 and 0.1)';
%          [Cat,Msg]=VO.GALEX.run_galex_sql(Q,[],'cell');
% Reliable: 
%------------------------------------------------------------------------------
TmpFile = tempname;

Dir       = Util.files.which_dir('run_galex_sql');
CasJobsDir = sprintf('%s%s..%sbin%sCasJobs%sGALEX%s',Dir,filesep,filesep,filesep,filesep,filesep);
CasJobsCL = sprintf('java-1.6.0 -jar casjobs.jar execute');
CasJobsCL = sprintf('java -jar /home/eran/matlab/fun/bin/CasJobs/GALEX/casjobs.jar execute');
CasJobsCL = sprintf('%scasjobs.jar execute',CasJobsDir);


if (nargin==1)
   OutputFileName = [];
   QueryType    = 'file';
   ConcatString = '';
elseif (nargin==2)
   QueryType    = 'file';
   ConcatString = '';
elseif (nargin==3)
   ConcatString = '';
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end


if (iscell(Query)==1)
   % if first argument is acell array set type to 'cell'
   QueryType = 'cell';
end

switch QueryType
 case 'file'
    %--- Read query from file ---
    FID = fopen(Query,'r');
    QueryLine = '';
    while (feof(FID)==0)
       Line       = fgetl(FID);
    
       %--- remove comment line started with # ---
       if (isempty(Line)==0)
          if (strcmp(Line(1),'#')),
             % comment line - ignore
          else
             QueryLine  = sprintf('%s %s',QueryLine,Line);
          end
       else
          %--- ignore empty lines ---
       end
    end
 case 'string'
    %--- Query is already in string format ---
    QueryLine = Query;
 case 'cell'
    % construct query from a cell array:

    % SELECT
    QueryLine = 'SELECT ';
    QueryLine = [QueryLine, fprintf_cell([],'%s, ',Query{1})];
    % remove last coma from string
    QueryLine = QueryLine(1:end-2);

    % FROM
    QueryLine = [QueryLine, ' FROM ', fprintf_cell([],'%s, ',Query{2})];
    % remove last coma from string
    QueryLine = QueryLine(1:end-2);

    % WHERE
    QueryLine = [QueryLine, ' WHERE ', Query{3}];
    % remove last coma from string

 otherwise
    error('Unknown QueryType Option');
end

%--- concatenate string to query ---
QueryLine = sprintf('%s %s',QueryLine,ConcatString);

if (isempty(OutputFileName)==1),
   OutputFileName = TmpFile;
   DeleteFile = 'y';
else
   DeleteFile = 'n';
end

%--- running SQL Query ---
%CL = sprintf('%s "%s" > %s',CasJobsCL,QueryLine,OutputFileName);
CL = sprintf('%s "%s"',CasJobsCL,QueryLine)
Status = 0;

PWD = pwd;
cd(CasJobsDir);
try
   [Status,Res] = system(CL);
catch
   Status = 1;   % failed
end
cd(PWD);


Msg = 0;
if (Status==1),
   % Failed in running python script
   Msg = 0;
else
   if (nargout>0),

      Str = str_duplicate('%f',length(Query{1})-1,'%f');
      Cell=textscan(Res,Str,'Delimiter',',','Headerlines',3);
      Cat = cell2mat(Cell);
      I = find(Cat==-999);
      Cat(I) = NaN;

      if (size(Cat,1)>0),
         Msg = 1;
      end
   end
end

