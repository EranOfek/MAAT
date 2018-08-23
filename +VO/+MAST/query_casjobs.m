function [Out,ColCell,Status,ResultOrig,Query]=query_casjobs(Query,varargin)
% Query MAST CasJobs service (requires casjobs.jar)
% Package: VO.MAST
% Description: Query the MAST CasJobs service (requires casjobs.jar)
%              For casjobs.jar instellation see:
%              http://mastweb.stsci.edu/mcasjobs/casjobscl.aspx
% Input  : - Query string or 3 elements query clause of the:
%            select, from, where parts of the query.
%            (e.g., {{'ra','dec'},{'photoobj'},''}).
%            Default is 'select top 10 ra, dec from photoobj'
%            If string start with '@' then assume the query is in a file.
%            If empty use default.
%          * Arbitrary number of ...,key,val,... pairs.
%            The following keywords are available:
%            'Table' - Query table name. Default is 'GALEX_GR6Plus7'.
%                      Other options include
%                      'PanSTARRS_DR1','HSCv2','GAIA_DR1','GALAEX_GR6Plus7',...
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
%            'SaveInTable' - Save catalog in table (rather than matrix).
%                      Default is false.
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
%            'WaitRetry' - Wait before retry [s]. Default is 1.
%            'Boolean2num' - convert boolean strings to numbers.
%                       Default is true.
%            'ShowResult' - Dipslay result string. Default is false.
%            'Verbose'- Verbose: true|false. Default is true.
% Output : - Output catalog
%          - Cell array of column names
%          - Result status.
%          - Result string.
%          - Query string.
% See also: query_casjobs_mydb
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Out,ColCell,Status]=VO.MAST.query_casjobs;
%          Out = VO.MAST.query_casjobs('select top 10 ra, dec from photoobj');
%          Out = VO.MAST.query_casjobs('select ra, dec from photoobj','boxcoo',[0 0.1 0 0.1].*pi./180,'StrRA','RA','StrDec','Dec');
%          Out = VO.MAST.query_casjobs('select top 10 raMean, decMean from ObjectThin','Table','PanSTARRS_DR1');
%          Out = VO.MAST.query_casjobs('select raMean, decMean from ObjectThin WHERE raMean>0 and raMean<0.1 and decMean>0 and decMean<0.1','Table','PanSTARRS_DR1');
%          % To get all table names
%          [~,~,~,Out] = VO.MAST.query_casjobs('SELECT table_name FROM information_schema.tables;')
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==0)
    Query = [];
end
if (isempty(Query))
    % set default query
    %Query = 'select * from information_schema.tables';
    Query = 'select top 10 ra, dec from photoobj';
end


DefV.Table                = 'GALEX_GR6Plus7';
DefV.CasJobsPath          = '~/matlab/bin/CasJobs/jar';
DefV.ConfigFilePath       = '~/matlab/bin/CasJobs/jar/CasJobs.config';
DefV.FormatString         = [];
DefV.OutType              = 'AstCat';
DefV.SaveInTable          = false;
DefV.SortBy               = [];
DefV.OnlyQuery            = false;
DefV.BoxCoo               = [];
DefV.StrRA                = [];
DefV.StrDec               = [];
DefV.WaitRetry            = 1;
DefV.Boolean2num          = true;
DefV.ShowResult           = false;
DefV.Verbose              = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% construct SQL query
Query = Util.sql.construct_query(Query,[],InPar.BoxCoo,InPar.StrRA,InPar.StrDec);

if (InPar.Verbose)
    fprintf('%s\n',Query);
end

%CasJobsPath    = '~/matlab/bin/CasJobs/jar';
%ConfigFilePath = '~/matlab/bin/CasJobs/jar/CasJobs.config';
%Query = 'select top 10 ra, dec from photoobj';

if (InPar.OnlyQuery)
    % Return only the query
    Out        = [];
    ColCell    = [];
    Status     = [];
    ResultOrig = [];
else
    %--- Normal execuation ---
    % execute the query
    
    PWD = pwd;

    cd(InPar.CasJobsPath);
    SysCommand = sprintf('java -jar %s/casjobs.jar execute -t "%s" "%s"',InPar.CasJobsPath,InPar.Table,Query);

    [Status,Result] = system(SysCommand);
    if (InPar.ShowResult)
        Result
    end
    
    % search for some known problems and try to recover:
    Problem1 = strfind(Result,'There is already an open D');
    Problem2 = strfind(Result,'was deadlocked on lock resources with another process');
    Problem3 = strfind(Result,'Object reference not set to an instance of an object');
    Problem4 = strfind(Result,'nested exception is:');
    if ~(isempty(Problem1) && isempty(Problem2) && isempty(Problem3) && isempty(Problem4))
        % some kind of communication problem try again (once)
        if (InPar.Verbose)
            fprintf('Some kind of communication problem try again (once)\n');
        end
        pause(InPar.WaitRetry);  % wait 1s before trying again
        [Status,Result] = system(SysCommand);
    end

    cd(PWD);
    ResultOrig = Result;

    Iend = strfind(Result,'Query complete!');
    if (isempty(Iend))
        %--------------------
        %--- query failed ---
        %--------------------
        if (InPar.Verbose)
            fprintf('Query failed\n');
            fprintf('Status = %d\n',Status);
            fprintf('%s\n',Result);
        end
        Status = -100;
        Out     = [];
        ColCell = {};
    
        if (~isempty(strfind(Result,'Query results exceed memory limit')))
            % Query exceed memory limit
            Status  = -101;
            if (InPar.Verbose)
                fprintf('Query results exceed memory limit\n');
            end
        end
    else
        %-----------------------
        %--- Query completed ---
        %-----------------------
        % check number of columns and their names
        RE       = regexp(Result,'\[\w+\]','match');    % get column names
        ColCell  = regexprep(RE,{'[',']'},'');          % remove brackets from column names
        Ncol     = numel(ColCell);                      % number of columns
        
        if (isempty(InPar.FormatString))
            InPar.FormatString = Util.string.str_duplicate('%f ',Ncol,'\n');
        end

        Result = Result(1:Iend-3);
        Istart = strfind(Result,'[');
        Result = Result(Istart(1):end);
        
        % deal with boolean columns
        if (InPar.Boolean2num)
            Result = regexprep(Result,'True','1');
            Result = regexprep(Result,'False','0');
        end

        C = textscan(Result,InPar.FormatString,'Delimiter',',','Headerlines',1);
        if (InPar.SaveInTable)
            Cat = table(C{:});
            Cat.Properties.VariableNames = ColCell;
        else
            Cat = cell2mat(C);
        end
        
        if (~isempty(Cat))
            switch lower(InPar.OutType)
                case 'astcat'
                    Out = AstCat;
                    CatField     = AstCat.CatField;
                    ColCellField = AstCat.ColCellField;
                    Out.(CatField) = Cat;
                    Out.(ColCellField) = ColCell;
                    Out = colcell2col(Out);
                case 'mat'
                    % do nothing
                    Out = Cat;
                otherwise
                    error('Unknown OutType option');
            end

            if (~isempty(InPar.SortBy))
                % sort output catalog
                Out = sortrows(Out,InPar.SortBy);
            end
        else
            Out = [];
        end

    end
end

