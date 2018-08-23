function SCat=astcat_search(CatName,RA,Dec,Radius,varargin)
% Search an astronomical catalog formatted in HDF5/HTM/zones format
% Example: VO.search.astcat_search('FIRST',1,1,0.01)


InvRAD = pi./180;
RAD    = 180./pi;


DefV.CooUnits             = 'rad';

DefV.DecZonesNameFormat   = '%s_GridHTMdecZones.hdf5';
DefV.DecZonesVarName      = '/V';
DefV.NdecZones            = 180;
DefV.GridHTMNameFormat    = '%s_GridHTM.hdf5';
DefV.GridHTMVarName       = '/V';
DefV.HTMNameFormat        = '%s_HTM_%07d.hdf5';
DefV.VarNameFormat        = '/HTM_%07d';
DefV.NvarInFile           = 1000;
DefV.Shape                = 'circ';  % 'none'|'circ'

if (~isempty(varargin))
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
else
    InPar = DefV;
end

ConvertCoo = convert.angular(InPar.CooUnits,'rad');
RA  = RA.*ConvertCoo;
Dec = Dec.*ConvertCoo;

LatStep = (180./InPar.NdecZones).*InvRAD;   % [radian]


DecZonesName = sprintf(InPar.DecZonesNameFormat,CatName);
GridHTMName  = sprintf(InPar.GridHTMNameFormat,CatName);

% Read HTM declination zone catalog
DecZone = h5read(DecZonesName,InPar.DecZonesVarName);


MinDec = max(Dec - Radius-LatStep,-pi./2);
MaxDec = min(Dec + Radius+LatStep,pi./2);

% DecZone is small
% direct search is faster than Util.find.bin_sear
F=find(MinDec<=DecZone(:,1) & MaxDec>=DecZone(:,1));

% tic;
% I1 = Util.find.bin_sear(DecZone(:,1),MinDec);
% I2 = Util.find.bin_sear(DecZone(:,1),MaxDec);
% toc

HTM_Level = 5;
HTM_Side  = 0.5.*pi./(2.^(HTM_Level - 1));
MinHTMind = min(DecZone(F,2));
MaxHTMind = max(DecZone(F,3));
NcolGridHTM = 16;
GridHTM = h5read(GridHTMName,InPar.GridHTMVarName,[MinHTMind, 1],[MaxHTMind-MinHTMind,NcolGridHTM]);
%tic;search_cat(GridHTM(:,2:3),RA,Dec); toc %,'SearchRad',3./RAD);toc
Ind = VO.search.search_sortedlat(GridHTM(:,2:3),RA,Dec,Radius+HTM_Side);
GridHTM   = GridHTM(Ind,:);
PolesLong = GridHTM(:,[11 13 15])';
PolesLat  = GridHTM(:,[12 14 16])';

Flag=celestial.htm.cone_in_polysphere(PolesLong,PolesLat,RA,Dec,Radius).' & GridHTM(:,4)>0.1;
HTMind = GridHTM(Flag,1);
NfinalCat = sum(GridHTM(Flag,4));
NcolCat   = 26;
Col.RA    = 1;
Col.Dec   = 2;
CatCooUnits  = 'rad';
ConvertCat   = convert.angular(CatCooUnits,'rad');

SCat = zeros(NfinalCat,NcolCat);
Nhtm = numel(HTMind);
K = 1;
for Ihtm=1:1:Nhtm
    FileInd = floor(HTMind(Ihtm)./InPar.NvarInFile).*InPar.NvarInFile;
    FileName = sprintf(InPar.HTMNameFormat,CatName,FileInd);
    VarName  = sprintf(InPar.VarNameFormat,HTMind(Ihtm));
    %Cat(Ihtm).Cat = h5read(FileName,VarName).';
    Dhtm = h5read(FileName,VarName);
    SizeD = size(Dhtm,1);
    SCat(K:K-1+SizeD,:) = Dhtm;
    K = K+SizeD;
end
if (Nhtm==0)
    SCat = [];
else
    switch lower(InPar.Shape)
        case 'none'
            % do nothing
        case 'circ'
            % within circle
            % sort and binary search
            
            % simple - search all
            D = celestial.coo.sphere_dist_fast(SCat(:,Col.RA).*ConvertCat,SCat(:,Col.Dec).*ConvertCat,RA,Dec);
            SCat = SCat(D<=Radius,:);
        otherwise
            error('Unknown Shape option');
    end
end

