function search_htmcat(CatName,Long,Lat,Radius,varargin)
% Search a local HTM/HDF5 catalog
% Description: Perform a cone search in an local catalog formatted in
%              the HTM (Hiraricical Triangular Mesh) grid and stired in
%              HDF5 format.
%              Such catalogs are constructed and described in the
%              Catalog.prep.build_htm_catalog function.
% Input  : - 
% Output : -


RAD = 180./pi;
ARCSEC_DEG = 3600;

CatName = 'FIRST';
Long    = pi;
Lat     = 0;
Radius  = 3600./(RAD.*3600);

UseZone = false;
MethodSearchGrid = 'find';
NvarInFile   = 1000;

% DefV.FilesDir             = './';
% DefV.ColCell              = [];
% DefV.ColUnits             = [];
% DefV.ColRA                = 1;
% DefV.ColDec               = 2;
% DefV.HTM_Level            = 9;
% DefV.HTM                  = [];
% DefV.LevList              = [];
% DefV.CatName              = '';
% DefV.DecZonesNameFormat   = '%s_GridHTMdecZones.hdf5';
% DefV.DecZonesVarName      = '/V';
% DefV.GridHTMNameFormat    = '%s_GridHTM.hdf5';
% DefV.GridHTMVarName       = '/V';
% DefV.HTMNameFormat        = '%s_HTM_%06d.hdf5';
% DefV.VarNameFormat        = '/HTM_%06d';
% DefV.NvarInFile           = 1000;
% DefV.NdecZones            = 180;
% DefV.CatRegion            = [0 360 -90 90]./RAD;  % [MinRA, MaxRA,MinDec,MaxDec] for data [radians]
% DefV.PrepCatGridHTM       = true;
% DefV.PrepCatGridHTMdecZones = true;
% DefV.SortHTM              = true;
% InPar = InArg.populate_keyval(DefV,varargin,mfilename);


DecZonesNameFormat   = '%s_GridHTMdecZones.hdf5';
DecZonesVarName      = '/V';
GridHTMNameFormat    = '%s_GridHTM.hdf5';
GridHTMVarName       = '/V';
HTMNameFormat        = '%s_HTM_%06d.hdf5';
VarNameFormat        = '/HTM_%06d';

RadiusRad = Radius./(RAD.*ARCSEC_DEG);
MinLat    = Lat-RadiusRad;
MaxLat    = Lat+RadiusRad;


GridHTMFile = sprintf(GridHTMNameFormat,CatName);

if (UseZone)
    %---------------------------------------------------
    % Search the object coordinates in the DecZone file
    %---------------------------------------------------
    DecZonesFile = sprintf(DecZonesNameFormat,CatName);
    %tic;Util.IO.loadh(DecZonesFile,DecZonesVarName);toc
    ZoneData = h5read(DecZonesFile,DecZonesVarName);
    Flag = MinLat>=ZoneData(1:end-1,1) & MaxLat<=ZoneData(2:end,1);
    MinHTMind = ZoneData(Flag,2);
    MaxHTMind = ZoneData(Flag,3);
    Count     = MaxHTMind - MinHTMind;

    % read the partial HTM grid file
    GridData    = h5read(GridHTMFile,GridHTMVarName,[MinHTMind 1],[Count 9]);
else
    % Read the entire HTM grid file
    GridData    = h5read(GridHTMFile,GridHTMVarName);
end

% Search the GridData file
% The columns in *_GridHTM.hdf5 file
ColInd  = 1;
ColRA   = 2;
ColDec  = 3;
ColNsrc = 4;

HTMcircumRadius = 10./RAD;
GridSearchRadius = HTMcircumRadius + RadiusRad;

switch lower(MethodSearchGrid)
    case 'find'
        Dist     = celestial.coo.sphere_dist_fast(GridData(:,ColRA),GridData(:,ColDec),Long,Lat);
        GridData = GridData(GridData(:,ColNsrc)>0 & Dist<=GridSearchRadius,:);
    case 'binary'
        I1 = bin_sear(GridData(:,ColDec),MinLat);
        I2 = bin_sear(GridData(:,ColDec),MaxLat);
        GridData = GridData(I1:I2,:);
        Dist     = celestial.coo.sphere_dist_fast(GridData(:,ColRA),GridData(:,ColDec),Long,Lat);
        GridData = GridData(GridData(:,ColNsrc)>0 & Dist<=GridSearchRadius,:);
    otherwise
        error('Unknown MethodSearchGrid option');
end

IndHTM  = GridData(:,ColInd);
IndFile = floor(IndHTM./NvarInFile);

%----------------------
% Search the HTM files 
%----------------------
Data  = [];
Nfile = size(GridData,1);
for Ifile=1:1:Nfile
    HTMFile = sprintf(HTMNameFormat,CatName,IndFile(Ifile));
    VarName = sprintf(VarNameFormat,IndHTM(Ifile));
    Data    = h5read(HTMFile,VarName);
end


Data