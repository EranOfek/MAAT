function [Out,ColCell,Status,ResultOrig,Query]=query_casjobs_recur(Query,varargin)
% Query MAST CasJobs service recursively for a box (requires casjobs.jar)
% Package: VO.MAST
% Description: Query the MAST CasJobs service recursively for a coordinates
%              in a box (requires casjobs.jar). For non recursive query use
%              VO.MAST.query_casjobs.
% Input  : - Query string or 3 elements query clause of the:
%            select, from, where parts of the query.
%            (e.g., {{'ra','dec'},{'photoobj'},''}).
%            Default is 'select top 10 ra, dec from photoobj'
%            If string start with '@' then assume the query is in a file.
%            If empty use default.
%          * Arbitrary number of ...,key,val,... pairs.
%            The following keywords are available:
%            'Table' - Query table name. Default is 'GALEX_GR6Plus7'.
%                      Other options include 'PanSTARRS_DR1','HSCv2',...
%            'CasJobsPath' - The CasJobs java file path
%                      The Java binary is available from:
%                      http://mastweb.stsci.edu/ps1casjobs/casjobscl.aspx
%                      Default is '~/matlab/bin/CasJobs/jar'.
%                      The configuration file named 'CasJobs.config'
%                      should be located in the same directory.
%            'FormatString' - Format string for query output.
%                      If empty, will attempt to use '%f %f...'.
%                      Default is [].
%            'OutType' - Output type 'AstCat'|'mat'.
%                      Default is 'AstCat'.
%            'SortBy' - Sort output by column index or name.
%                      Default is [].
%            'OnlyQuery' - A flag indicating if to return only the query.
%                      If truen, then the query string will be returned
%                      without execuation. Default is false.
%            'BoxCoo' - Optional [MinRA, MaxRA, MinDec, MaxDec] to query.
%                       Default is empty [radians].
%                       If not empty must also provide 'StrRA' and
%                       'StrDec'. 
%            'StrRA'  - RA string name. Default is empty.
%            'StrDec' - Dec string name. Default is empty.
%            'NumSubBox' - Recursive coordinate division. Default is 3.
%            'Verbose'- Verbose: true|false. Default is true.
% Output : - Output catalog
%          - Cell array of column names
%          - Result status.
%          - Result string.
%          - Query string.
% See also: query_casjobs_mydb
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Out = VO.MAST.query_casjobs_recur('select ra, dec from photoobj','boxcoo',[0 0.1 0 0.1].*pi./180,'StrRA','ra','StrDec','dec');
%          Out = VO.MAST.query_casjobs_recur('select raMean, decMean from ObjectThin','Table','PanSTARRS_DR1','boxcoo',[0 0.5 0 0.5].*pi./180,'StrRA','raMean','StrDec','decMean');
% Reliable: 2
%--------------------------------------------------------------------------

CatField = AstCat.CatField;

if (nargin==0)
    Query = [];
end
if (isempty(Query))
    % set default query
    %Query = 'select * from information_schema.tables';
    Query = 'select top 10 ra, dec from photoobj';
end


DefV.Table                = 'GALEX_GR6Plus7';
DefV.CasJobsPath          = '~/matlab/bin/CasJobs/GALEX';
DefV.ConfigFilePath       = '~/matlab/bin/CasJobs/GALEX/CasJobs.config';
DefV.FormatString         = [];
DefV.OutType              = 'AstCat';
DefV.SortBy               = [];
DefV.OnlyQuery            = false;
DefV.BoxCoo               = [];
DefV.StrRA                = [];
DefV.StrDec               = [];
DefV.Recursion            = true;
DefV.NumSubBox            = 3;
DefV.Verbose              = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% Execute the query non recursively
[Out,ColCell,Status,ResultOrig]=VO.MAST.query_casjobs(Query,varargin{:});

if (Status==-101)
    % Query failed due to "exceed memory limit"
    % attempt to decrease the box search size and run recursively
    
    DeltaRA  = InPar.BoxCoo(2) - InPar.BoxCoo(1);
    DeltaDec = InPar.BoxCoo(4) - InPar.BoxCoo(3);

    fprintf('Recursive box query\n');

    K = 0;
    B = zeros(InPar.NumSubBox.^2,4);
    for Isize=0:1:InPar.NumSubBox-1
        for Jsize=0:1:InPar.NumSubBox-1
            K = K + 1;
            B(K,:) = [Isize, Isize+1, Jsize, Jsize+1];
        end
    end
    ListEdge = B./InPar.NumSubBox;
    Nedge = size(ListEdge,1);

    switch lower(InPar.OutType)
        case 'mat'
            Out = [];
        case 'astcat'
            Out = AstCat;
        otherwise
            error('Unknown OutMat option');
    end
    
    for Iedge=1:1:Nedge

        fprintf('Recursive box query for sub box %d\n',Iedge);

        MinRA1   = InPar.BoxCoo(1)  + ListEdge(Iedge,1).*DeltaRA;
        MaxRA1   = InPar.BoxCoo(1)  + ListEdge(Iedge,2).*DeltaRA;
        MinDec1  = InPar.BoxCoo(3)  + ListEdge(Iedge,3).*DeltaDec;
        MaxDec1  = InPar.BoxCoo(3)  + ListEdge(Iedge,4).*DeltaDec;

        % call query_casjobs recursively
        [Cat,ColCell,Status,ResultOrig]=VO.MAST.query_casjobs_recur(Query,varargin{:},'BoxCoo',[MinRA1, MaxRA1, MinDec1, MaxDec1]);
        
        % Merge results
        if (isnumeric(Cat))
            Out = [Out;Cat];
        else
            if (isempty(Out.(CatField)))
                Out = Cat;
            else
                Out = merge([Out;Cat]);
            end
        end
    end
end

    
if (~isempty(Out))

    if (~isempty(InPar.SortBy))
        % sort output catalog
        Out = sortrows(Out,InPar.SortBy);
    end
end