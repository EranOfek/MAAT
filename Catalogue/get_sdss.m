function [Cat,ColCell,Col]=get_sdss(RA,Dec,SearchRad,Shape,CatType)
%--------------------------------------------------------------------------
% get_sdss function                                              Catalogue
% Description: Search by coordinates local copy of the SDSS-DR10 source
%              catalog.
%              OBSOLETE: Use VO.search.get_sdss instead.
% Input  : - 
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    May 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cat,ColCell,Col]=get_sdss(180./RAD,0./RAD,0.01./RAD);
% Reliable: 
%--------------------------------------------------------------------------
InvRAD = pi./180;

Def.SearchRad = 0.5.*InvRAD;   % radians
Def.Shape     = 'circ';
Def.CatType   = 'struct';
if (nargin==2),
    SearchRad = Def.SearchRad;
    Shape     = Def.Shape;
    CatType   = Def.CatType;
elseif (nargin==3),
    Shape     = Def.Shape;
    CatType   = Def.CatType;
elseif (nargin==4),
    CatType   = Def.CatType;
elseif (nargin==5),
    % do nothing
else
    error('Illegal number of input arguments');
end

% V and Verr units are [deg/day]
ColCell = {'ra','dec','type','flags','modelMag_u','modelMagErr_u','modelMag_g','modelMagErr_g','modelMag_r','modelMagErr_r','modelMag_i','modelMagErr_i','modelMag_z','modelMagErr_z','V','Verr'};
Ncol    = numel(ColCell);
Col     = cell2struct(num2cell(1:1:Ncol),ColCell,2);


%DefV. = 
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

RA  = convertdms(RA,'gH','R');
Dec = convertdms(Dec,'gD','r');

Nlev         = 10;
TriSide      = 180./(2.^Nlev);              % estimated size of the triangle side
TriCM        = TriSide.*(sqrt(2)./2).*1.01;  
SearchRadHTM = SearchRad + TriCM.*InvRAD;

%SDSS_htmCat = load_check('SDSS_htmCat.mat');
SDSS_htmCat = loadh('SDSS_htmCat.hd5','SDSS_htmCat');
CatH = search_cat(SDSS_htmCat,RA,Dec,'CooType','sphere','SearchMethod','binms','SearchRad',SearchRadHTM);
IndHTM = SDSS_htmCat(CatH.IndCat,3);
Nch  = numel(IndHTM);
Cat = zeros(0,Ncol);

for Ich=1:1:Nch,
     FN    = sprintf('SDSS_htmS%07d.hd5',floor(IndHTM(Ich)./1000).*1000);
     SubFN = sprintf('SDSS_htm%07d',IndHTM(Ich));
     Cat   = [Cat;loadh(FN,SubFN)];
end

D = sphere_dist(Cat(:,1),Cat(:,2),RA,Dec);
Cat = Cat(D<SearchRad,:);

switch lower(CatType)
    case 'struct'
        Tmp = Cat;
        clear Cat;
        Cat.Cat     = Tmp;
        Cat.ColCell = ColCell;
        Cat.Col     = Col;
    otherwise
        % do nothing
end

