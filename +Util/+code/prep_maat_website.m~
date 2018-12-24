function prep_maat_website(varargin)
% Prepare the Matlab Astronomy & AStrophysics Toolbox website
% Package: Util.code
% Description: Prepare the Matlab Astronomy & AStrophysics Toolbox website
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Util.code.prep_maat_website
% Reliable: 1
%--------------------------------------------------------------------------

DefV.FunPath             = '~/matlab/fun/';
DefV.DataPath            = '~/matlab/data/';
DefV.WebPath             = '/raid/eran/public_html/matlab/';
DefV.DocPath             = '/raid/eran/public_html/matlab/doc/';
DefV.WebFunPath          = '/raid/eran/public_html/matlab/matlab/fun/';
DefV.ManualPath          = '~/matlab/fun/+manual/';
DefV.MainTemplateFile    = 'MainTemplate.html';
DefV.FunTableName        = 'FunTable.html';
DefV.FunListName         = 'FunList.html';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

PWD = pwd;

% read main html template file
cd(InPar.WebPath);
MainTemplateStr = Util.files.file2str(sprintf('%s',InPar.MainTemplateFile),'str');


% prepare function list
cd(InPar.FunPath);
[List,TableCell]=Util.code.prep_function_list;
% convert date to YYYY-mm-DD
DD=datestr(TableCell(:,4),'YYYY-mm-DD');
for I=1:1:numel(List)
    TableCell{I,4} = DD(I,:);
end

% generate HTML table of functions
cd(InPar.WebPath);
Header = {'Function Name','Package/Class','IsStatic','Last update','Brief description'};
www.html_table(InPar.FunTableName,'TableCell',TableCell,'TableStatment','class="sortable"','TableHead',Header);

% add sortable javascript
FileStr = Util.files.file2str(InPar.FunTableName,'str');
Comment = sprintf('<b>Table is sortable - click on column to sort.</b><br> IsStatic indicate if the function is a static class.<br> Number of functions: %d.<br> Last update: %s.<br>',numel(List),date);
FileStr = sprintf('<script src="sorttable.js"></script>\n\n %s\n<br>\n %s',Comment,FileStr);
FID = fopen(InPar.FunTableName,'w');
fprintf(FID,'%s',FileStr);
fclose(FID);


FileStr = regexprep(MainTemplateStr,'PUT_HERE_HTML_FILENAME',InPar.FunTableName);
FID = fopen(InPar.FunListName,'w');
fprintf(FID,'%s',FileStr);
fclose(FID);






% convert all mlx files in manual directory to html files
cd (InPar.ManualPath);
% FilesMlx = dir('*.mlx');
% Nf = numel(FilesMlx);
% for If=1:1:Nf
%     publish(FilesMlx(If).name,'html');
% end

% copy all html files in +manual directory to doc/ dir
FilesHtml = dir('*.html');
Nf = numel(FilesHtml);
for If=1:1:Nf
    Dest = sprintf('%s%s_Content.html',InPar.DocPath,FilesHtml(If).name(1:end-5));
    copyfile(FilesHtml(If).name,Dest);
end

% generate a main web page for each html document file
cd(InPar.DocPath)
%FilesCont = dir('*_Content.html');
FilesContC = regexprep({FilesHtml.name},'.html','_Content.html');
Nc = numel(FilesContC);
[FilesCont(1:Nc).name] = deal(FilesContC{:});

Nf = numel(FilesCont);
for If=1:1:Nf
    FileName = sprintf('%s.html',FilesCont(If).name(1:end-13));
    FileName
    FileNameContent = sprintf('%s_Content.html',FilesHtml(If).name(1:end-5));
    FileStr = regexprep(MainTemplateStr,'PUT_HERE_HTML_FILENAME',['./doc/',FileNameContent]);
    FID = fopen(FileName,'w');
    fprintf(FID,'%s',FileStr);
    fclose(FID);
    
    
end

% add html href target
for If=1:1:Nf
    FileName = sprintf('%s.html',FilesCont(If).name(1:end-13));
    FileNameContent = sprintf('%s_Content.html',FilesHtml(If).name(1:end-5));
    system(sprintf('sed ''s/<a href/<a target="_parent" href/g'' < %s > %s1',FileNameContent,FileNameContent));
    system(sprintf('mv %s1 %s',FileNameContent,FileNameContent));
end

% copy documentation to index.html
system('cp ./documentation.html index.html');

% copy main_Content.html
cd(InPar.WebPath);
system('mv ./doc/main.html index.html');
%system('mv ./doc/main_Content.html .');




% rsync function directory
cd(InPar.FunPath)
system(sprintf('rsync -avx --exclude=".*" ./ %s',InPar.WebFunPath));


% prepare tar files
Itar = 0;
cd(InPar.FunPath)
cd('../')
% tar file of fun dir
Itar = Itar + 1;
TarFilesData(Itar).FileName   = 'Fun.tar';
TarFilesData(Itar).FileNameGZ = sprintf('%s.gz',TarFilesData(Itar).FileName);
system(sprintf('tar -cvf %s fun/*',TarFilesData(Itar).FileName)); 
system(sprintf('gzip %s',TarFilesData(Itar).FileName));
DirFun  = dir(sprintf('%s',TarFilesData(Itar).FileNameGZ));
TarFilesData(Itar).FileSize = DirFun.bytes./1024./1000;
system(sprintf('mv %s %s',TarFilesData(Itar).FileNameGZ,InPar.WebPath));
% tar file for data/+cats
Itar = Itar + 1;
TarFilesData(Itar).FileName   = 'data_cats.tar';
TarFilesData(Itar).FileNameGZ = sprintf('%s.gz',TarFilesData(Itar).FileName);
system(sprintf('tar -cvf %s data/+cats/*',TarFilesData(Itar).FileName)); 
system(sprintf('gzip %s',TarFilesData(Itar).FileName));
DirFun  = dir(sprintf('%s',TarFilesData(Itar).FileNameGZ));
TarFilesData(Itar).FileSize = DirFun.bytes./1024./1000;
system(sprintf('mv %s %s',TarFilesData(Itar).FileNameGZ,InPar.WebPath));
% tar file for ...


% startup.m
system(sprintf('cp startup.m %s',InPar.WebPath));
cd('fun/+Util/+code/');
system(sprintf('cp install_maat.m %s',InPar.WebPath));

% copy to kesem
cd(InPar.WebPath);
input('press enter to continue','s');
system('rsync -avx ./ eofek@kesem:/ph1users/eofek/public_html/matlab/');


% return to original directory
cd(PWD);