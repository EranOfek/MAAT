function build_htm4cat(CatName,Level)
%--------------------------------------------------------------------------
% build_htm4cat function                                         Catalogue
% Description: Given an astronomical catalog (in a structure format,
%              e.g., 'FIRST.mat') construct a directory containing
%              catalogs per HTM (Hierarchical Triangular Mesh) region,
%              and an index catalog of all HTM regions.
% Input  : - Catalog name for which to construct an HTM version
%            (e.g., 'NVSS.mat').
%          - HTM level of the catalog.
% Output : null
% See also: htm_build.m, search_htm_coocat.m, get_apass.m
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: build_htm4cat('FIRST.mat',6)
% Reliable: 2
%--------------------------------------------------------------------------

DirSuffix = '_htm';
[HTM,LevList]=celestial.htm.htm_build(Level);

CatAll = Util.IO.load2(CatName);
CatBaseName = regexp(CatName,'.mat','split');
CatBaseName = CatBaseName{1};
DirName = sprintf('%s%s%s',CatBaseName,DirSuffix,filesep);


CatDir  = Util.files.which_dir(CatName);
PWD = pwd;
cd(CatDir);
mkdir(DirName);

LevelList = LevList(end);
Nptr      = numel(LevelList.ptr);
MasterCat.Cat = zeros(Nptr,4);
for Il=1:1:Nptr,
    Ptr = LevList(Level).ptr(Il);
    Flag = celestial.htm.in_polysphere(CatAll.Cat(:,[CatAll.Col.RA, CatAll.Col.Dec]),HTM(Ptr).cosd);
    
    Cat      = CatAll.Cat(Flag,:);
    Cat      = sortrows(Cat,CatAll.Col.Dec);
    %Cat.Col      = CatAll.Col;
    %Cat.ColCell  = CatAll.ColCell;
    %Cat.ColUnits = CatAll.ColUnits;
    %Cat.SortedBy = 'Dec';
    %Cat.SortedByCol = CatAll.Col.Dec;
    if (size(Cat,1)>0),
        save(sprintf('%s%s%s%05d.mat',DirName,CatBaseName,DirSuffix,Il),'Cat');
    end
    
    %MeanCD = mean(HTM(Ptr).cosd);
    %[MeanRA,MeanDec] = celestial.coo.cosined2coo(MeanCD(1),MeanCD(2),MeanCD(3));
    %[MeanRA,MeanDec] = tri_equidist_center(HTM(Ptr).coo(:,1).', HTM(Ptr).coo(:,2).');
    MeanRA = mean(HTM(Ptr).coo(:,1).');
    MeanDec = mean(HTM(Ptr).coo(:,2).');
    MasterCat.Cat(Il,:) = [MeanRA, MeanDec, Il, size(Cat,1)];
end

Cat = MasterCat;
Cat.Cat = sortrows(Cat.Cat,2);
Cat.ColCell  = {'RA','Dec','Ptr','Nsrc'};
Cat.Col      = cell2struct(num2cell(1:1:length(Cat.ColCell)),Cat.ColCell,2);
Cat.ColUnits = {'rad','rad','',''};
Cat.SortedBy = 'Dec';
Cat.SortedByCol = 2;
Cat.Description = sprintf('Catalog of %s HTM pointers',CatBaseName);
Cat.isHTM       = true;
% CatHTM field
Cat.CatHTM.Col         = CatAll.Col;
Cat.CatHTM.ColCell     = CatAll.ColCell;
Cat.CatHTM.ColUnits    = CatAll.ColUnits;
Cat.CatHTM.SortedBy    = 'Dec';
Cat.CatHTM.SortedByCol = CatAll.Col.Dec;
Cat.FileBaseName       = sprintf('%s%s%%05d.mat',CatBaseName,DirSuffix);

save(sprintf('%s%s.mat',CatBaseName,DirSuffix),'Cat');
cd(PWD);

