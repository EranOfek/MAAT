function Cat=search_mhtm_cat(CatName,RA,Dec,SearchRadius,varargin)
% Search master HTM catalog 
% Package: Catalog.search
% Description: A Master HTM catalog contains a list of HTM regions and the
%              HDF file name containing the HTM region. This function can
%              search for objects near coordinates using such master HTM
%              catalogs.
% Input  : - Master HTM catalog name in HDF5 format.
%          - J2000.0 R.A. 
%          - J2000.0 Dec.
%          - Search radius [arcsec default]. Default is 60".
%          * Arbitrary number of ...,key,val,.. arguments.
%            The following keywords are available:
%            'SearcRadUnits' - Search radius units. See convert.angular for
%                          options. Default is 'arcsec'.
%            'BaseName'  - HTM files base name.
%                          Default is 'GAIA_DR1_HTM%06d.hdf5'.
%            'VarBaseName'- HDF5 variable name base name.
%                          Default is 'HTM%06d'.
%            'HTMinFile' - Number of HTM in file. Default is 100.
%            'ColCellFile' - Mat file containing ColCell of catalog column
%                          names. Default is 'GAIA_DR1_HTMcolcell.mat'.
%            'SortBy'    - Sort catalog by column name or index.
%                          If empty then do not sort.
%                          Default is [].
%            'ColRA'     - RA column name in catalog. Default is 'RA'.
%            'ColDec'     - Dec column name in catalog. Default is 'Dec'.
% Example: Cat=Catalog.search.search_mhtm_cat('GAIA_DR1_HTMindex.hdf5',1,1,3600);
% Reliable: 2

RAD = 180./pi;
ARCSEC_DEG = 3600;

CatField  = AstCat.CatField;
ColHTM.Index = 3;

DefV.SearcRadUnits        = 'arcsec';
DefV.BaseName             = 'GAIA_DR1_HTM%06d.hdf5';
DefV.VarBaseName          = 'HTM%06d';
DefV.HTMinFile            = 100;
DefV.ColCellFile          = 'GAIA_DR1_HTMcolcell.mat';
DefV.SortBy               = [];
DefV.ColRA                = 'RA';
DefV.ColDec               = 'Dec';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% convert search radius to arcsec
SearchRadius = SearchRadius.*convert.angular(InPar.SearcRadUnits,'arcsec');

% convert RA/Dec to rad
RA  = celestial.coo.convertdms(RA,'gH','r');
Dec = celestial.coo.convertdms(Dec,'gD','R');

% Read the Master HTM catalog
%MasterHTM = AstCat.loadh2astcat(CatName);

Tmp = Util.IO.loadh(CatName);
FN = fieldnames(Tmp);
MasterHTM.(CatField) = Tmp.(FN{1});

% Number of HTMs
Nhtm      = size(MasterHTM.(CatField),1);
% Estimate HTM search radius for bounding circle
Level  = floor(log(Nhtm)./log(4));
RadHTM = sqrt(2).*90./(2.^Level)./RAD;   % [rad]

% Search Master HTM catalog
% assume catalog is sorted by Dec
Res = search_cat(MasterHTM.(CatField),RA,Dec,'SearchRad',SearchRadius./(RAD.*ARCSEC_DEG)+RadHTM);

% Index of selected entries in CatHTM
SelectedInd = [Res.IndCat];
SelectedInd = MasterHTM.(CatField)(SelectedInd,ColHTM.Index);


% load all individual HDF5 catalogs
Nsel = numel(SelectedInd);
if (Nsel==0)
    % nothing found
    Cat = AstCat;
else
    
    for Isel=1:1:Nsel
        % catalog HTM index
        IndHTM  = SelectedInd(Isel);
        IndFile = floor(double(IndHTM)./InPar.HTMinFile).*InPar.HTMinFile;
        % construct catalog name
        Name = sprintf(InPar.BaseName,IndFile);
        VarName = sprintf(InPar.VarBaseName,IndHTM);

        % load catalog
        C(Isel).Cat = Util.IO.loadh(Name,VarName).';

    end
    % merge catalogs
    Cat         = AstCat;
    Cat.Cat     = [C.(CatField)].';
    Cat.ColCell = Util.IO.load2(InPar.ColCellFile);
    Cat         = colcell2col(Cat);
    
    % sort
    if (~isempty(InPar.SortBy))
        Cat = sortrows(Cat,InPar.SortBy);
    end
    
    % remove sources outside search region
    ColRA  = colname2ind(Cat,InPar.ColRA);
    ColDec = colname2ind(Cat,InPar.ColDec);
    CatField = AstCat.CatField;
    
    %Dist = sphere_dist(Cat,RA,Dec,{'RA','Dec'},'rad');
    Dist = celestial.coo.sphere_dist_fast(Cat.(CatField)(:,ColRA),Cat.(CatField)(:,ColDec),RA,Dec);

    %Dist.Dist = sphere_dist(Cat.Cat(:,1),Cat.Cat(:,2),RA,Dec);
    %Flag = Dist.Dist<SearchRadius./(RAD.*ARCSEC_DEG);
    Flag = Dist<SearchRadius./(RAD.*ARCSEC_DEG);
    
    Cat.(CatField)  = Cat.(CatField)(Flag,:);

    % clean catalog

end



