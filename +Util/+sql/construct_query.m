function Query=construct_query(Query,Into,BoxCoo,StrRA,StrDec)
% Construct an SQL query from SELECT, FROM, WHERE clauses.
% Package: Util.sql
% Description: Construct an SQL query from SELECT, FROM, WHERE clauses.
% Input  : - Query string or 3 elements query clause of the:
%            select, from, where parts of the query.
%            (e.g., {{'ra','dec'},{'photoobj'},''}).
%            Default is 'select top 10 ra, dec from photoobj'
%            If string start with '@' then assume the query is in a file.
%            If empty use default.
%          - If not empty then an "into %s" string will be added.
%            Default is empty.
%          - Optional Box coordinates [MinRA, MaxRA, MinDec, MaxDec] in
%            radians.
%          - RA string, must be provided if box coordinates are given.
%          - Dec string, must be provided if box coordinates are given.
% Output : - Query string.
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Query=Util.sql.construct_query({{'RA','Dec'},{'someTable'},'RA>0'},'mydb.test')
%          Query=Util.sql.construct_query({{'RA','Dec'},{'someTable'},''},'mydb.test',[0 1 0 1].*pi./180,'RA','Dec')
% Reliable: 2


RAD = 180./pi;

if (nargin<2)
    Into = [];
    BoxCoo = [];
else
    if (nargin<3)
        BoxCoo = [];
    end
end

AddAnd = false;
if (iscell(Query))
    % construct query from a cell array:

    if iscell(Query{1})
        % SELECT
        QueryLine = 'SELECT ';
        QueryLine = [QueryLine, Util.IO.fprintf_cell([],'%s, ',Query{1})];
        % remove last coma from string
        QueryLine = QueryLine(1:end-2);
    else
        QueryLine = Query{1};
    end

    if iscell(Query{2})
        % FROM
        QueryLine = [QueryLine, ' FROM ', Util.IO.fprintf_cell([],'%s, ',Query{2})];
        % remove last coma from string
        QueryLine = QueryLine(1:end-2);
    else
        QueryLine = sprintf('%s %s',QueryLine,Query{2});
    end

    % WHERE
    if (~isempty(Query{3}))
        AddAnd = true;
    end
    QueryLine = [QueryLine, ' WHERE ', Query{3}];
    % remove last coma from string
        
    
    Query = QueryLine;
else
    if (strcmp(Query(1),'@'))
        % assume query is in a file
        QueryFile = Query(2:end);
        Query     = Util.files.file2str(QueryFile,'str');
    end
end

% add box coordinates
if (~isempty(BoxCoo))
    Iwh = strfind(lower(Query),'where');
    if (isempty(Iwh))
        Query = sprintf('%s WHERE',Query');
    end
    if (AddAnd)
        Query = sprintf('%s and',Query');
    end
    BoxCoo = BoxCoo.*RAD;
    %Query  = sprintf('%s %s between %f and %f and %s between %f and %f',Query,StrRA,BoxCoo(1),BoxCoo(2),StrDec,BoxCoo(3),BoxCoo(4));
    Query  = sprintf('%s %s>%f and %s<%f and %s>%f and %s<%f',Query,StrRA,BoxCoo(1),StrRA,BoxCoo(2),StrDec,BoxCoo(3),StrDec,BoxCoo(4));
end


% Into option:
% save the file in the local MAST disk
if (~isempty(Into))
    Query = sprintf('%s into %s',Query,Into);
end
