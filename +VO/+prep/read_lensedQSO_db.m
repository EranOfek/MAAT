function T=read_lensedQSO_db(varargin)
% Read garvitationaly lensed quasars database
% Package: VO.prep
% Description: Read garvitationaly lensed quasars database from:
%              https://www.ast.cam.ac.uk/ioa/research/lensedquasars/index.html
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'URL' - default is 'https://www.ast.cam.ac.uk/ioa/research/lensedquasars/index.html'
% Output : - AstCat object containing the list of lensed quasars.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: T=VO.prep.read_lensedQSO_db;
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;


DefV.URL                  = 'https://www.ast.cam.ac.uk/ioa/research/lensedquasars/index.html';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


[Cell,Header] = www.parse_html_table_old(InPar.URL);

N = numel(Cell{1});
Tmp = regexp(Cell{2},'https://www.ast.cam.ac.uk/ioa/research/lensedquasars/indiv/(?<url>.+\**\.html)\''>(?<name>[\w\+\-\*]+)</a>','names');


URL  = cell(N,1);
Name = cell(N,1);
for I=1:1:N
    URL{I}  = sprintf('https://www.ast.cam.ac.uk/ioa/research/lensedquasars/indiv/%s',Tmp{I}.url);
    Name{I} = Tmp{I}.name;
end

RA  = str2double(Cell{3});
Dec = str2double(Cell{4});

W1W2  = str2double(Cell{6});
zQSO  = str2double(Cell{7});
zLens = str2double(Cell{8});
Nim   = str2double(Cell{9});
GaiaDen = str2double(Cell{10});
G     = str2double(Cell{15});
Disc  = Cell{end};

T = AstCat;
T.Cat = table(RA,Dec,Name,W1W2,zQSO,zLens,Nim,GaiaDen,G,Disc,URL);
T.Cat.RA  = T.Cat.RA./RAD;
T.Cat.Dec = T.Cat.Dec./RAD;

T.Cat.Properties.VariableNames = {'RA','Dec','Name','W1W2','zQSO','zLens','Nim','GaiaDen','G','Disc','URL'};
T.ColCell = T.Cat.Properties.VariableNames;
T.ColUnits = {'rad','rad','','mag','','','','','mag','',''};
T = colcell2col(T);
T.Source = InPar.URL;
