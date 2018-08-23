function [Res,Dist,PA]=sdss_dr2_search(RA,Dec,ObjType,SearchSize,SearchShape);
%---------------------------------------------------------------------------------------------
% sdss_dr2_search function            Load and search the local copy compact version of 
%                                   the SDSS DR2 cataloge.
%                                   The Catalogue contains some basics info (only) on 
%                                   each star or galaxy in the SDSS DR2 catalog.
% Input  : - R.A. if scalar than given in radinas, if three element vector than [H M S].
%          - Dec. if scalar than given in radinas, if four element vector than [Sign D M S].
%          - ObjType: {'Gal' | 'Star'}.
%          - Search radius/half-width in arcsec.
%          - Shape of search region, {'circle','box'}, default is 'box'.
% Output : - Catalog of objects within the search radius.
%            If empty than no sources has been found.
%            Columns:
%   ra, dec, primTarget,
%   psfMag_u, extinction_u, psfMagErr_u,
%   psfMag_g, extinction_g, psfMagErr_g,
%   psfMag_r, extinction_r, psfMagErr_r,
%   psfMag_i, extinction_i, psfMagErr_i,
%   psfMag_z, extinction_z, psfMagErr_z
%          - Vector of respective distances [radians] between search coordinates and object.
%          - Vector of respective P.A. [radians] between search coordinates and object.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    July 2004     
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------------------------
RAD               = 180./pi;
ARCSEC_DEG        = 3600;
STARCAT_ROOT_NAME = 'SDSS_DR2_Star_';
GALCAT_ROOT_NAME  = 'SDSS_DR2_Gal_';
EXTCAT_ROOT_NAME  = '.mat';
CAT_STEP_SIZE     = 0.05;

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
    Res            = zeros(0,18);
 case 'Gal'
    FILE_ROOT_NAME = GALCAT_ROOT_NAME;
    Res            = zeros(0,18);
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
