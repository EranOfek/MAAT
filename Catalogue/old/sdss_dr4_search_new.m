function [Res2,Dist2,PA2]=sdss_dr4_search_temp_test(RA,Dec,ObjType,SearchSize,SearchShape,TileList);
%---------------------------------------------------------------------------------------------
% sdss_dr4_search function            Load and search the local copy compact version of 
%                                   the SDSS DR3 cataloge.
%                                   The Catalogue contains some basics info (only) on 
%                                   each star or galaxy in the SDSS DR3 catalog.
% Input  : - R.A. if scalar than given in radinas, if three element vector than [H M S].
%          - Dec. if scalar than given in radinas, if four element vector than [Sign D M S].
%          - ObjType: {'Gal' | 'Star' | 'Sec'}.
%          - Search radius/half-width in arcsec.
%          - Shape of search region:
%            'circle', default.
%            'box' - NOT WORKIG YET!
% Output : - Catalog of objects within the search radius.
%            If empty than no sources has been found.
%            Columns for Star catalog:
%            (1)     - RA   [radians]
%            (2)     - Dec  [radians]
%            (3)     - type
%            (4)     - flags
%            (5)     - status
%            (6)     - primTarget
%            (7)     - probPSF
%            (8-12)  - run, rerun, camcol, field, obj,
%            (13-15) - psfMag_u, extinction_u, psfMagErr_u,
%            (16-18) - psfMag_g, extinction_g, psfMagErr_g,
%            (19-21) - psfMag_r, extinction_r, psfMagErr_r,
%            (22-24) - psfMag_i, extinction_i, psfMagErr_i,
%            (25-27) - psfMag_z, extinction_z, psfMagErr_z,
%            (28-32) - texture_u, texture_g, texture_r, texture_i, texture_z,
%            (33-37) - lnLStar_u, lnLStar_g, lnLStar_r, lnLStar_i, lnLStar_z,
%            (38-41) - mE1_i, ME2_i, mRrCc_i, mCr4_i,
%            (42-44) - isoA_i, isoB_i, isoPhi_i,
%            (45-48) - rowv, rowvErr, colv, colvErr
%
%            Columns for Gal catalog:
%            (1)     - RA   [radians]
%            (2)     - Dec  [radians]
%            (3)     - type
%            (4)     - flags
%            (5)     - status
%            (6)     - primTarget
%            (7)     - probPSF
%            (8-12)  - run, rerun, camcol, field, obj,
%            (13-15) - modelMag_u, extinction_u, modelMagErr_u,
%            (16-18) - modelMag_g, extinction_g, modelMagErr_g,
%            (19-21) - modelMag_r, extinction_r, modelMagErr_r,
%            (22-24) - modelMag_i, extinction_i, modelMagErr_i,
%            (25-27) - modelMag_z, extinction_z, modelMagErr_z,
%            (28-32) - fracDeV_u, fracDeV_g, fracDeV_r, fracDeV_i, fracDeV_z,
%            (33-38) - isoA_r, isoAErr_r, isoB_r, isoBErr_r, isoPhi_r, isoPhiErr_r,
%            (39-44) - deVRad_r, deVRadErr_r, deVAB_r, deVABErr_r, deVPhi_r, deVPhiErr_r
%            (45-50) - deVRad_i, deVRadErr_i, deVAB_i, deVABErr_i, deVPhi_i, deVPhiErr_i
%            (51-56) - expRad_r, expRadErr_r, expAB_r, expABErr_r, expPhi_r, expPhiErr_r
%            (57-62) - expRad_i, expRadErr_i, expAB_i, expABErr_i, expPhi_i, expPhiErr_i
%
%            Columns for Sec catalog:
%            (1)     - RA   [radians]
%            (2)     - Dec  [radians]
%            (3)     - type
%            (4)     - flags
%            (5)     - status
%            (6)     - primTarget
%            (7)     - probPSF
%            (8-12)  - run, rerun, camcol, field, obj,
%            (13-15) - modelMag_u, extinction_u, modelMagErr_u,
%            (16-18) - modelMag_g, extinction_g, modelMagErr_g,
%            (19-21) - modelMag_r, extinction_r, modelMagErr_r,
%            (22-24) - modelMag_i, extinction_i, modelMagErr_i,
%            (25-27) - modelMag_z, extinction_z, modelMagErr_z,
%            (28-31) - rowv, rowvErr, colv, colvErr
%
%          - Vector of respective distances [radians] between search coordinates and object.
%          - Vector of respective P.A. [radians] between search coordinates and object.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    September 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------------------------
RAD               = 180./pi;
ARCSEC_DEG        = 3600;
STARCAT_ROOT_NAME = 'stars/SDSS_DR4_PrimStar_';
GALCAT_ROOT_NAME  = 'galaxies/SDSS_DR4_PrimGal_';
SECCAT_ROOT_NAME  = 'SDSS_DR4_Sec_';
EXTCAT_ROOT_NAME  = '.mat';
CAT_STEP_SIZE     = 0.05;

N_RA              = 360;
N_DEC             = 180;
N_SUB_TILES       = [5 5];

DEC_SIGN_STR      = ['n';'p'];

NcolStar          = 48;     % number of columns in Star cat
NcolGal           = 62;     % number of columns in Gal cat
NcolSec           = 31;     % number of columns in Sec cat

ColRA             = 1;
ColDec            = 2;

if (nargin==4),
   SearchShape = 'box';
elseif (nargin==6),
   % do nothing
else
   error('Ilegal number of input arguments');
end

if (length(RA)==1),
   % radians
elseif (length(RA)==3),
   RA  = convertdms(RA,'H','r');
else
   error('Illegal number of elements in RA');
end

if (length(Dec)==1),
   % radians
elseif (length(Dec)==4),
   Dec  = convertdms(Dec,'D','R');
else
   error('Illegal number of elements in RA');
end

switch ObjType
 case 'Gal'
    RootName = GALCAT_ROOT_NAME;
    NcolOut  = NcolGal;
    Type=sprintf('PrimGal');
 case 'Star'
    RootName = STARCAT_ROOT_NAME;
    NcolOut  = NcolStar;
    Type=sprintf('PrimStar');
 case 'Sec'
    RootName = SECCAT_ROOT_NAME;
    NcolOut  = NcolSec;
    Type=sprintf('Sec');
 otherwise
    error('Unknown ObjType Option');
end


% search size in radians
SearchSizeRad = SearchSize./(ARCSEC_DEG.*RAD); 

switch SearchShape
 case 'box'
    error('box search not available');    
    SearchRadFactor = sqrt(2);
 case 'circle'
    SearchRadFactor = 1;
 otherwise
    error('Unknown SearchShape Option');
end

% listing the relavent tiles
RA_max_search=RA+SearchSize.*(1./(ARCSEC_DEG.*RAD));
RA_min_search=RA-SearchSize.*(1./(ARCSEC_DEG.*RAD));
Dec_min_search=Dec-SearchSize.*(1./(ARCSEC_DEG.*RAD));
Dec_max_search=Dec+SearchSize.*(1./(ARCSEC_DEG.*RAD));

TilesInds=find((TileList(:,3)<=RA_max_search)&(TileList(:,4)>=RA_min_search)&(TileList(:,5)<=Dec_max_search)&(TileList(:,6)>=Dec_min_search))

Ntiles       = size(TilesInds,1)

N_SubTiles=5;

%--- Load all Tiles ---
Res = zeros(0,NcolOut);
for I=1:1:Ntiles,
    
   CenterRA   = TileList(TilesInds(I),1);
   CenterDec  = TileList(TilesInds(I),2);
   MinRA      = TileList(TilesInds(I),3);
   MaxRA      = TileList(TilesInds(I),4);
   MinDec     = TileList(TilesInds(I),5);
   MaxDec     = TileList(TilesInds(I),6);

   for IsubRA=1:1:N_SubTiles,
            for IsubDec=1:1:N_SubTiles,
                
               CenterSubRA  = MinRA  + (0.5+IsubRA -1).*(MaxRA - MinRA)./N_SubTiles;
               MinSubRA     = MinRA  + (0.0+IsubRA -1).*(MaxRA - MinRA)./N_SubTiles;
               MaxSubRA     = MinRA  + (1.0+IsubRA -1).*(MaxRA - MinRA)./N_SubTiles;
   
               CenterSubDec = MinDec + (0.5+IsubDec-1).*(MaxDec-MinDec)./N_SubTiles;
               MinSubDec    = MinDec + (0.0+IsubDec-1).*(MaxDec-MinDec)./N_SubTiles;
               MaxSubDec    = MinDec + (1.0+IsubDec-1).*(MaxDec-MinDec)./N_SubTiles;
               
               if ((MinSubRA<=RA_max_search)&(MaxSubRA>=RA_min_search)&(MinSubDec<=Dec_max_search)&(MaxSubDec>=Dec_min_search))

                   %--- Prep Tile variable ---
                   if (CenterDec>=0),
                       SignStr = 'p';
                   else
                       SignStr = 'n';
                   end
                   VarName = sprintf('/home/castor/assafh/tmp/SDSS/dr4/%s%05d_%s%04d',RootName,floor(CenterSubRA.*RAD.*100),SignStr,floor(abs(CenterSubDec).*RAD.*100));
                   TileFileName = sprintf('%s.mat',VarName);


                   try
                       load(TileFileName);
                   catch
                       % File doesn't exist
                       % continue
                       '?'
                       Cat = zeros(0,NcolOut);
                   end

                   Res = [Res; Cat];

               end
               
            end
            
   end   
   
end
%pause
%RA
%Dec
%size(Res)
%ColRA
%Res(1451,1);
%pause
[Dist,PA] = sphere_dist(RA,Dec,Res(:,ColRA),Res(:,ColDec));
ind=find(Dist<=SearchSize.*(1./(ARCSEC_DEG.*RAD)));
Res2=Res(ind,:);
Dist2=Dist(ind,:);
PA2=PA(ind,:);
Res(:,1:2);
%TileFileName














