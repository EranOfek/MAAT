function [Cat,Msg,FirstLine,Col,ColCell]=run_sdss_sql(Query,varargin) 
%--------------------------------------------------------------------------
% run_sdss_sql function                                               sdss
% Description: Run SQL query on SDSS database and retrieve the results
%              into an array.
% Input  : - File name containing the SQL query, or alternatively,
%            a string containing the full query,  or alternatively,
%            this can be a 3 elements cell array, in which the 
%            first cell contains : a cell array of attributes in the
%                                  SELECT clause.
%            second cell contains: a cell array of tables/views in the
%                                  FROM clause.
%            third cell contains : a string to appear in the WHERE clause. 
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'OutputFileName' - File name. If empty then no output file.
%                               Default is empty.
%            'QueryType' - One of the following query type:
%                          'file' - First argument is file name
%                                   containing the query.
%                                   Default (unless query is a cell).
%                          'string' - First argument is the full query
%                                   string.
%                          'cell'   - First argument is a cell array of
%                                   attributes.
%            'ConcatString' - Optional string to concatenate to query
%                       (e.g., this could be a variable part of the query),
%                       default is empty string (i.e., '').
%            'ServerURL' - Default is 'http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
% Output : - Optional catalog, containing the search result.
%          - Return code: 1 - objects found; 0- No objects have been found
%          - First line in resulted file.
%          - Structure containing a list of all column indices.
%          - Cell array of column names.
% See also: wget_sdss.m
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Static file query:
%          Query='/scr2/eran/bin/sdssQA/CasJOB/sql/get_1004.sql';
%          [Cat,Msg,FirstLine]=run_sdss_sql(Query);
%          % dynamic file query:
%          Query='/scr2/eran/matlab/sdss/bin/sql/get_obj_field.sql';
%          Con='run=3813 and rerun=41 and camcol=6 and field=60';
%          [Cat,Msg,FirstLine]=run_sdss_sql(Query,'QueryType','file','ConcatString',Con);
%          % string query:
%          QueryStr = 'select top 2 ra,dec from star';
%          [Cat,Msg,FirstLine]=run_sdss_sql(QueryStr,'QueryType','string');
%          % Cell query:
%          Q{1} = {'ra','dec','psfMag_g'};
%          Q{2} = {'Star'};
%          Q{3} = '(ra between 0 and 1) and (dec between 0 and 1)';
%          [Cat,Msg,FirstLine]=run_sdss_sql(Q);
% Reliable: 1
%--------------------------------------------------------------------------
import Util.IO.*

DefV.OutputFileName = [];
DefV.QueryType      = 'file';
DefV.ConcatString   = '';
DefV.OutFormats     = 'csv'; %['csv','xml','html']
%DefV.ServerURL      = 'http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx';
%DefV.ServerURL      = 'http://skyserver.sdss.org/dr12/en/tools/search/x_sql.aspx';
%DefV.ServerURL      = 'http://skyserver.sdss.org/dr12/en/tools/search/sql.aspx';
%DefV.ServerURL      = 'http://skyserver.sdss.org/SkyserverWS/dr12/SearchTools/SqlSearch';
DefV.ServerURL      = 'http://skyserver.sdss.org/dr12/en/tools/search/x_results.aspx';

% astro_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
% public_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (iscell(Query)==1)
   % if first argument is acell array set type to 'cell'
   InPar.QueryType = 'cell';
end

switch InPar.QueryType
 case 'file'
    %--- Read query from file ---
    FID = fopen(Query,'r');
    QueryLine = '';
    while (feof(FID)==0)
       Line       = fgetl(FID);
    
       %--- remove comment line started with # ---
       if (isempty(Line)==0)
          if (strcmp(Line(1),'#'))
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
StrQuery = sprintf('%s %s',QueryLine,InPar.ConcatString);

%StrQuery = 'SELECT ra, dec, psfMag_g FROM Star WHERE (ra between 0 and 0.1) and (dec between 0 and 0.1)';
%Output    = urlread(InPar.ServerURL,'post',{'cmd',StrQuery,'format',InPar.OutFormats});
Output    = urlread(InPar.ServerURL,'post',...
    {'cmd',StrQuery,'format',InPar.OutFormats,'searchtool','SQL','TaskName','Skyserver.Search.SQL','syntax','NoSyntax','ReturnHtml','false'});
Sp        = regexp(Output,'\n','split');
FirstLine = Sp{2};
ColCell    = regexp(FirstLine,',','split');
ColCell   = regexprep(ColCell,',','_');
Format    = '';
for I=1:1:length(ColCell)
    Format = sprintf('%s %%f',Format);
end
Cat = textscan(Output,Format,'Delimiter',',','headerlines',2);
Cat = cell2mat(Cat);

Msg = '';

if (nargout>3)
    Col = cell2struct(num2cell((1:1:length(ColCell))),ColCell,2);
end
