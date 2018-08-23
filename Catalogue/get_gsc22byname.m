function Cat=get_gsc22byname(Name);
%-----------------------------------------------------------------------------
% get_gsc22byname function                                          Catalogue
% Description: Get GSC-2.2 catalog from VizieR web server, by GSC2 name.
% Input  : - GSC2 name.
% Output : - Structure containing catalog:
%            .Data - 20 column cell array with the retrieved objects.
%            .Header - 20 Column cell array with the description of each column.
%            .Units  - 20 Column cell array with the units of each column.
%            .Col    - Structure with the column number for each field.
%            Comments:
%                 Class: 0=star 1=Galaxy, 2=Blend 3=Non-star 4=unclassified 5=Defect
% Tested : Matlab 7.3
%     By : Eran O. Ofek            November 2008
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: C=get_gsc22byname('N1202223119');
%-----------------------------------------------------------------------------
DirProg    = Util.files.which_dir('get_gsc22');   % local directory
DirViz     = sprintf('%s%s%s%s%s','bin',filesep,'vizquery',filesep,'cdsclient-2.87');
ProgName   = 'findgsc2.2';
MaxN       = 10000;
TmpFile    = 'tmp.gsc2.2';



Prog     = sprintf('%s%s..%s%s%s%s',DirProg,filesep,filesep,DirViz,filesep,ProgName);



ExecStr =  sprintf('%s -g %s > %s',Prog,Name,TmpFile);

delete(TmpFile);
system(ExecStr);

TmpFile1 = sprintf('%s_1',TmpFile);
delete(TmpFile1);
system(sprintf('sed "s#,# #g" %s | sed "s#;# #g" > %s',TmpFile,TmpFile1));
%system(sprintf('fgrep -v aclient %s > %s',TmpFile,TmpFile1));
delete(TmpFile);

FID = fopen(TmpFile1,'r');
Scan = textscan(FID,'%s %f %f %f    %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Headerlines',1,'CommentStyle','#');

Cat.Data = Scan;
Cat.Header{1} = 'GSC2-ID';
Cat.Header{2} = 'RA (J2000)';
Cat.Header{3} = 'Dec (J2000)';
Cat.Header{4} = 'Epoch';
Cat.Header{5} = 'eRA';
Cat.Header{6} = 'eDec';
Cat.Header{7} = 'Fmag';
Cat.Header{8} = 'eFmag';
Cat.Header{9} = 'Jmag';
Cat.Header{10} = 'eJmag';
Cat.Header{11} = 'Vmag';
Cat.Header{12} = 'eVmag';
Cat.Header{13} = 'Nmag';
Cat.Header{14} = 'eNmag';
Cat.Header{15} = 'Class';
Cat.Header{16} = 'Axis';
Cat.Header{17} = 'Ecc';
Cat.Header{18} = 'PA';
Cat.Header{19} = 'Status';
Cat.Header{20} = 'Dist';

Cat.Units{1} = '';
Cat.Units{2} = 'deg';
Cat.Units{3} = 'deg';
Cat.Units{4} = 'year';
Cat.Units{5} = 'arcsec';
Cat.Units{6} = 'arcsec';
Cat.Units{7} = 'mag';
Cat.Units{8} = 'mag';
Cat.Units{9} = 'mag';
Cat.Units{10} = 'mag';
Cat.Units{11} = 'mag';
Cat.Units{12} = 'mag';
Cat.Units{13} = 'mag';
Cat.Units{14} = 'mag';
Cat.Units{15} = '';
Cat.Units{16} = 'pix';
Cat.Units{17} = '';
Cat.Units{18} = 'deg';
Cat.Units{19} = '';
Cat.Units{20} = 'arcsec';

Cat.Col.ID     = 1;
Cat.Col.RA     = 2;
Cat.Col.Dec    = 3;
Cat.Col.Epoch  = 4;
Cat.Col.eRA    = 5;
Cat.Col.eDec   = 6;
Cat.Col.Fmag   = 7;
Cat.Col.eFmag  = 8;
Cat.Col.Jmag   = 9;
Cat.Col.eJmag  = 10;
Cat.Col.Vmag   = 11;
Cat.Col.eVmag  = 12;
Cat.Col.Nmag   = 13;
Cat.Col.eNmag  = 14;
Cat.Col.Class  = 15;
Cat.Col.Axis   = 16;
Cat.Col.Ecc    = 17;
Cat.Col.PA     = 18;
Cat.Col.Status = 19;
Cat.Col.Dist   = 20;
