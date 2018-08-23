function [Lines,DistL,PAL]=icat_search(Catalog,RA,Dec,Radius);
%---------------------------------------------------------------------------
% icat_search function                                            Catalogue
% Description: Search for an object in an astronomical catalog, where the
%              catalog is specified by its name.
% Input  : - Catalog name - supported catalogue:
%            'FIRST'   - FIRST
%            'NVSS'    - NVSS
%            'ROSAT'   - ROSAT faint source catalog
%            'ROSATb'  - ROSAT bright source catalog
%          - Single RA (radians, [HMS] or sexagesimal string),
%            or string of object name to resolve in SIMBAD - in this
%            case the Dec field must be empty (i.e., []).
%          - Single Dec (radians, [HMS] or sexagesimal string).
%            If empty then assume object name to resolve is given instead
%            of RA.
%          - Search radius in arcsec.
% Output : - The selected lines from the catalogue.
%          - Distance [rad] from search coordinates.
%          - Position angle [rad] from search coordinates.
% See Also: cat_search.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                 Feb 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
RAD = 180./pi;

if (isempty(Dec)==1),
   [RA,Dec] = get_simbad_coo(RA);
end

RA  = convertdms(RA,'gH','r');
Dec = convertdms(Dec,'gD','R');

switch Catalog
 case 'FIRST'
    load first_030411.mat
    Cat = first_030411;
    clear first_030411;
    ColRA  = 1;
    ColDec = 2;
    CatSortCol = 'Dec';
 case 'NVSS'
    load NVSS.mat;
    Cat = NVSS;
    clear NVSS;
    ColRA  = 1;
    ColDec = 2;
    CatSortCol = 'Dec';
 case 'ROSAT'
    load rosat_faint_xs.mat;
    Cat = rosat_faint_xs;
    clear rosat_faint_xs.mat;
    Cat(:,1:2) = Cat(:,1:2)./RAD;
    ColRA  = 1;
    ColDec = 2;
    CatSortCol = 'Dec';
 case 'ROSATb'
    load rosat_bsc_1rxs.mat;
    Cat = rosat_bsc_1rxs;
    clear rosat_bsc_1rxs;
    Cat(:,1:2) = Cat(:,1:2)./RAD;
    ColRA  = 1;
    ColDec = 2;
    CatSortCol = 'Dec';
 otherwise
    error('Unsupported catalog name');
end

[Ind,DistL,PAL]=cat_search(Cat,[ColRA, ColDec],[RA, Dec], Radius./(3600.*RAD), 'circle', CatSortCol);

Lines = Cat(Ind,:);
