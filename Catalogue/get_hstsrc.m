function Cat=get_hstsrc(RA,Dec,Radius)
%--------------------------------------------------------------------------
% get_hstsrc function                                            Catalogue
% Description: Search a retrieve the HST source catalog around given
%              coordinate.
% Input  : - J2000.0 R.A. (radians, sexagesimal string or [H M S]).
%          - J2000.0 Dec. (radians, sexagesimal string or [Sign D M S]).
%          - Search radius [radians].
% Output : - Catalog of sources found in HST images.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Mar 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=get_hstsrc(148.8883./RAD,69.0653./RAD,0.02./RAD)
% Reliable: 2
%--------------------------------------------------------------------------

CatName_HSTsrc = 'HSTsrc.mat';
Path           = Util.files.which_dir(CatName_HSTsrc);
SubDir         = 'HSTsrc_ImID';
CatFileName    = 'HSTsrc_ImID%07d.mat';

%DefV. = 
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

HSTsrc = Util.IO.load2(CatName_HSTsrc);

Res = search_cat(HSTsrc.Cat(:,[HSTsrc.Col.RA, HSTsrc.Col.Dec]),RA,Dec,'SearchRad',Radius,'SearchMethod','binms');

PathCat = sprintf('%s%s%s%s',Path,filesep,SubDir,filesep);
Ncat    = numel(Res.IndCat);
AllCat  = zeros(0,13);
for Icat=1:1:Ncat
    Cat = Util.IO.load2(sprintf('%s%s',PathCat,sprintf(CatFileName,HSTsrc.Cat(Res.IndCat(Icat),HSTsrc.Col.ImageID) )));
    AllCat = [AllCat;Cat.Cat];
end
D = celestial.coo.sphere_dist(AllCat(:,Cat.Col.RA),AllCat(:,Cat.Col.Dec),RA,Dec);
Cat.Cat = AllCat(D<Radius,:);
Cat.Cat = sortrows(Cat.Cat,Cat.Col.Dec);

