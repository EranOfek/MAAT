function [OutFile,Status,Result]=query_casjobs_mydb(Query,varargin)
% Query MAST CasJobs service into MAST mydb (requires casjobs.jar)
% Package: VO.MAST
% Description: Query the MAST CasJobs service (requires casjobs.jar),
%              but store the output in mydb.
% Input  : - Query string or 3 elements query clause of the:
%            select, from, where parts of the query.
%            (e.g., {{'ra','dec'},{'photoobj'},''}).
%            Default is 'select top 10 ra, dec from photoobj'
%            If string start with '@' then assume the query is in a file.
%            If empty use default.
%          * Arbitrary number of ...,key,val,... pairs.
%            The following keywords are available:
%            'Table' - Query table name. Default is 'PanSTARRS_DR1'.
%                      Other options include 'GALEX_GR6Plus7','HSCv2',...
%            'Into'  - String indicating into which DB name to save the
%                      results. Default is 'test' which translate to.
%            'CasJobsPath' - The CasJobs java file path
%                      The Java binary is available from:
%                      http://mastweb.stsci.edu/ps1casjobs/casjobscl.aspx
%                      Default is '~/matlab/bin/CasJobs/jar'.
%                      The configuration file named 'CasJobs.config'
%                      should be located in the same directory.
%            'Extract' - A flag indicating if to extract the results from
%                      the MAST DB archive. Default is true.
%            'BaseURL' - MAST archive URL.
%                      Default is
%                      'http://mastweb.stsci.edu/CasOutPut/CSV/'.
%            'UserName' - MAST archive user name.
%                      Default is 'EranOfek'.
%            'Drop'   - Drop the table from the MAST archive after
%                      retrival. Default is true.
%            'Verbose'- Verbose: true|false. Default is true.
% Output : - The output file name in which the query results are stored.
%          - The execuation status.
%          - The execuation result string.
%            Note that the actual result is stored in the MAST mydb
%            database.
% See also: query_casjobs
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [F,St,Res]=VO.MAST.query_casjobs_mydb('select top 10 ra, dec from photoobj');
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


DefV.Table                = 'PanSTARRS_DR1';
DefV.Into                 = 'test';
DefV.CasJobsPath          = '~/matlab/bin/CasJobs/GALEX';
DefV.ConfigFilePath       = '~/matlab/bin/CasJobs/GALEX/CasJobs.config';
DefV.Extract              = true;
DefV.BaseURL              = 'http://mastweb.stsci.edu/CasOutPut/CSV/';
DefV.UserName             = 'EranOfek';
DefV.Drop                 = false;
DefV.Verbose              = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

MydbInto  = sprintf('mydb.%s',InPar.Into);
Query = Util.sql.construct_query(Query,[]); %MydbInto);


if (InPar.Verbose)
    fprintf('%s\n',Query);
end


PWD = pwd;

cd(InPar.CasJobsPath);
SysCommand = sprintf('java -jar %s/casjobs.jar execute -t "%s" "%s"',InPar.CasJobsPath,InPar.Table,Query);

[Status,Result] = system(SysCommand);

% search for some known problems and try to recover:
Problem1 = strfind(Result,'There is already an open D');
Problem2 = strfind(Result,'was deadlocked on lock resources with another process');
Problem3 = strfind(Result,'Object reference not set to an instance of an object');
if ~(isempty(Problem1) && isempty(Problem2) && isempty(Problem3))
    % some kind of communication problem try again (once)
    if (InPar.Verbose)
        fprintf('Some kind of communication problem try again (once)\n');
    end
    [Status,Result] = system(SysCommand);
end

OutFile = [];
if (InPar.Extract)
    % Example: casjobs extract -b row1 -F -type csv -d
    SysCommand = sprintf('java -jar %s/casjobs.jar extract -b %s -F -type csv -d',InPar.CasJobsPath,InPar.Into);
    [Status1,Result1] = system(SysCommand);

    OutFile = sprintf('%s_%s',InPar.Into,InPar.UserName);
    URL     = sprintf('%s%s',InPar.BaseURL,OutFile);
    
    cd(PWD);
    www.pwget(URL);
    cd(InPar.CasJobsPath);

    if (InPar.Drop)
        % Example: casjobs execute -t "PANSTARRS_DR1" -n "drop query" "drop table row1"
        SysCommand = sprintf('java -jar %s/casjobs.jar execute -t "%s" -n "drop query" "drop table %s"',InPar.CasJobsPath,InPar.Table,InPar.Into);
        [Status2,Result2] = system(SysCommand);
    end
end
 
cd(PWD);


