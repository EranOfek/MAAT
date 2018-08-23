function [Res,Dist,PA]=sdss_dr3_search(RA,Dec,ObjType,SearchSize,SearchShape);
%---------------------------------------------------------------------------------------------
% sdss_dr3_search function            Load and search the local copy compact version of 
%                                   the SDSS DR3 cataloge.
%                                   The Catalogue contains some basics info (only) on 
%                                   each star or galaxy in the SDSS DR3 catalog.
% Input  : - R.A. if scalar than given in radinas, if three element vector than [H M S].
%          - Dec. if scalar than given in radinas, if four element vector than [Sign D M S].
%          - ObjType: {'Gal' | 'Star'}.
%          - Search radius/half-width in arcsec.
%          - Shape of search region, {'circle','box'}, default is 'box'.
% Output : - Catalog of objects within the search radius.
%            If empty than no sources has been found.
%            Columns:
%            (1)     - ra [deg]
%            (2)     - dec [deg]
%            (3)     - primTarget,
%            (4-6)   - psfMag_u, extinction_u, psfMagErr_u,
%            (7-9)   - psfMag_g, extinction_g, psfMagErr_g,
%            (10-12) - psfMag_r, extinction_r, psfMagErr_r,
%            (13-15) - psfMag_i, extinction_i, psfMagErr_i,
%            (16-18) - psfMag_z, extinction_z, psfMagErr_z,
%            (19-23) - Run, ReRun, CamCol, Field, Object
%               For galaxy type also return:
%            (24-28) - fracDeV_u, fracDeV_g, fracDeV_r, fracDeV_i, fracDeV_z,
%            (29-34) - isoA_r, isoAErr_r, isoB_r, isoBErr_r, isoPhi_r, isoPhiErr_r,
%            (35-40) - deVRad_r, deVRadErr_r, deVAB_r, deVABErr_r, deVPhi_r, deVPhiErr_r,
%            (41-46) - deVRad_i, deVRadErr_i, deVAB_i, deVABErr_i, deVPhi_i, deVPhiErr_i,
%            (47-52) - expRad_r, expRadErr_r, expAB_r, expABErr_r, expPhi_r, expPhiErr_r,
%            (53-58) - expRad_i, expRadErr_i, expAB_i, expABErr_i, expPhi_i, expPhiErr_i
%
%          - Vector of respective distances [radians] between search coordinates and object.
%          - Vector of respective P.A. [radians] between search coordinates and object.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    July 2004     
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------------------------
RAD               = 180./pi;
ARCSEC_DEG        = 3600;
STARCAT_ROOT_NAME = 'SDSS_DR3_Star_';
GALCAT_ROOT_NAME  = 'SDSS_DR3_Gal_';
EXTCAT_ROOT_NAME  = '.mat';
CAT_STEP_SIZE     = 0.05;
NcolStar          = 23;
NcolGal           = 58;

if (nargin==4),
   SearchShape = 'box';
elseif (nargin==5),
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

% search size in radians
SearchSizeRad = SearchSize./(ARCSEC_DEG.*RAD); 


RA_deg   = RA.*RAD;
Dec_deg  = Dec.*RAD;

DecLower = Dec - SearchSizeRad;
DecUpper = Dec + SearchSizeRad;

DecLower_Region = floor(DecLower.*RAD./CAT_STEP_SIZE).*CAT_STEP_SIZE;
DecUpper_Region = floor(DecUpper.*RAD./CAT_STEP_SIZE).*CAT_STEP_SIZE;

ListRegionsDec = [DecLower_Region:0.05:DecUpper_Region].';

switch ObjType
 case 'Star'
    FILE_ROOT_NAME = STARCAT_ROOT_NAME;
    Res            = zeros(0,NcolStar);
 case 'Gal'
    FILE_ROOT_NAME = GALCAT_ROOT_NAME;
    Res            = zeros(0,NcolGal);
 otherwise
    error('Unknown ObjType option');
end

Dist = zeros(0,1);
PA   = zeros(0,1);

%--- search for objects ---
for Iregion=1:1:length(ListRegionsDec),
   if (ListRegionsDec(Iregion)<0),
      Hem = 's';
   else
      Hem = 'n';
   end
   VarName   = sprintf('%s%s%04d',FILE_ROOT_NAME, Hem, abs(floor(ListRegionsDec(Iregion).*100)));
   FileName  = sprintf('%s%s',VarName, EXTCAT_ROOT_NAME);

   if (exist(FileName,'file')==0),
      %--- region not found
   else
      load(FileName);
      eval(sprintf('Cat = %s;',VarName));
      eval(sprintf('clear %s;',VarName));
      %--- Cat is sorted by RA ---
      [Lines,DistL,PAL] = cat_search(Cat(:,1:2)./RAD,[1 2],[RA,Dec],SearchSizeRad,SearchShape,'RA');

      Res   = [Res; Cat(Lines,:)];
      Dist  = [Dist; DistL];
      PA    = [PA; PAL];
   end      
end
