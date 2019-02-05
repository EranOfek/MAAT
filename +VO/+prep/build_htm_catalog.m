function Nsrc=build_htm_catalog(Cat,varargin)
% Build an HTM catalog in HDF5 format for fast queries
% Package: @VO.prep
% Description: Given a catalog, prepare the catalog in HDF5 format sorted
%              into HTM cells. This program is used for the construction of
%              fast access large catalogs. You can search such catalogs
%              using the VO.search.htmcat_search.
% Input  : - Catalog, in array or AstCat format, with Long/Lat in radians. 
%          * Arbitrary number of pairs of ...,key,val,... parameters.
%            The following keywords are available:
%            'Nsrc' - a matrix of [IndHTM Nsrc], where IndHTM is the HTM
%                     index, and Nsrc is the number of sources in HTM.
%            'SaveCol' - Index of columns to save. If empty, then save all
%                     columns. Default is empty.
%            'DecRange' - Dec range in between to save HTMs (radians).
%                     Default is [-pi./2 pi./2]
%            'ColCell'  - Cell array of column names. Default is {}.
%            'ColUnits' - Cell array of column units. Default is {}.
%            'HTM_Level'- HTM level. Default is 7.
%            'ColRA'    - RA column.
%            'ColDec'   - Dec column.
%            'CatName'  - Catalog base name.
%            'NfilesInHDFf' - Number of HTM catalogs in each HDF file.
%                         Default is 100.
%            'IndStep'  - Step for index file. Default is 30.
%            'SaveInd'  - Save index HDF file. Default is true.
% Output : null
% Example: [HTM,LevList] = celestial.htm.htm_build(5);
%          F=Util.IO.load2('FIRST.mat')
%          VO.prep.build_htm_catalog(F,'CatName','FIRST','HTM_Level',6);
% Reliable: 2

RAD = 180./pi;


DefV.Nsrc                 = [];   % Nsrc matrix [IndHTM, Nsrc]
DefV.SaveCol              = [];   % index of columns to save []->all
DefV.DecRange             = [-pi./2 pi./2];
DefV.RARange              = [0 2.*pi];
DefV.ColCell              = {};
DefV.ColUnits             = {};
DefV.HTM_Level            = 7;
DefV.ColRA                = 1;
DefV.ColDec               = 2;
DefV.CatName              = 'FIRST';
DefV.NfilesInHDF          = 100;
DefV.IndStep              = 30;
DefV.SaveInd              = true;
DefV.HTM                  = [];
DefV.LevelHTM             = [];

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (AstCat.isastcat(Cat))
    InPar.ColCell  = Cat.ColCell;
    InPar.ColUnits = Cat.ColUnits;
    Cat            = Cat.Cat;
end


RadiusHTM = (sqrt(2).*90./(2.^(InPar.HTM_Level - 1)))./RAD;
Radius    = 0.00001./(RAD.*3600);

% select columns
if (~isempty(InPar.SaveCol))
    Cat = Cat(:,InPar.SaveCol);
end

% sort cat by declination
Cat = sortrows(Cat,InPar.ColDec);

% cut catalog
I1 = Util.find.bin_sear(Cat(:,InPar.ColDec),InPar.DecRange(1)-RadiusHTM);
I2 = Util.find.bin_sear(Cat(:,InPar.ColDec),InPar.DecRange(2)+RadiusHTM);
CatC    = Cat(I1:I2,:);

% build HTM
if (~isempty(InPar.HTM) && ~isempty(InPar.LevelHTM))
    HTM      = InPar.HTM;
    LevelHTM = InPar.LevelHTM;
else
    [HTM,LevelHTM] = celestial.htm.htm_build(InPar.HTM_Level);
end

ListIndexHTM   = LevelHTM(InPar.HTM_Level).ptr;
Nhtm           = numel(ListIndexHTM);

if (isempty(InPar.Nsrc))
    Nsrc = nan(Nhtm,2);
else
    Nsrc = InPar.Nsrc;
end

for Ihtm=1:1:Nhtm
    Ihtm
    % check if HTM mean Dec is in dec range
    IndHTM = ListIndexHTM(Ihtm);
    MeanRA  = mean(HTM(IndHTM).coo(:,1));
    MeanDec = mean(HTM(IndHTM).coo(:,2));
    
    if (MeanDec>=InPar.DecRange(1) && MeanDec<InPar.DecRange(2) && MeanRA>=InPar.RARange(1) && MeanRA<InPar.RARange(2))
        % HTM in dec/ra range
        
        
        D = celestial.coo.sphere_dist_fast(CatC(:,InPar.ColRA),CatC(:,InPar.ColDec),MeanRA,MeanDec);
        FlagD = D<(1.5.*LevelHTM(end).side);
        CatCC = CatC(FlagD,:);
        
        Flag = celestial.htm.in_polysphere(CatCC(:,[InPar.ColRA, InPar.ColDec]),HTM(IndHTM).coo,2);
        
        %Flag = celestial.htm.cone_in_polysphere(HTM(IndHTM).PolesCoo(:,1),HTM(IndHTM).PolesCoo(:,2),CatC(:,InPar.ColRA), CatC(:,InPar.ColDec),Radius);
        Nsrc(Ihtm,:) = [IndHTM, sum(Flag)];
        
        if (Nsrc(Ihtm,2)>0)
            [FileName,DataName]=HDF5.get_file_var_from_htmid(InPar.CatName,IndHTM,InPar.NfilesInHDF);
            HDF5.save_cat(FileName,DataName,CatCC(Flag,:),InPar.ColDec,InPar.IndStep);
        end
        
    end
    
end

% save HTM index file
if (InPar.SaveInd)
    IndFileName = sprintf('%s_htm.hdf5',InPar.CatName);
    delete(IndFileName);
    % Nsrc=HDF5.get_nsrc(CatName);
    HDF5.save_htm_ind(HTM,IndFileName,[],{},Nsrc)

    HDF5.save_cat_colcell(InPar.CatName,InPar.ColCell,InPar.ColUnits);
end
