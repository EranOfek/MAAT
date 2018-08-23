function [Out,QS]=wget_hsc_sources(RA,Dec,varargin)
% Query sources in the HST source catalog tables
% Package: VO.MAST
% Description: Query sources in the HST catalog tables. By default the
%              query is done on a predefined columns in the SExtractor
%              catalogs of the WFPC2, ACS, and WFPC3.
% Input  : - J2000.0 RA [radians, sexagesimal string or h,m,s]
%          - J2000.0 Dec [radians, sexagesimal string or sign,d,m,s]
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'HalfSize'  - Search box half size. Default is 10 arcsec.
%            'HalfSizeUnits' - Search half size units. Default is 'arcsec'.
%            'SearchTable'   - Tables in which to search for sources.
%                          Default is: {'Catalog_WFPC2_SourceExtractor','Catalog_ACS_SourceExtractor','Catalog_WFC3_SourceExtractor'};
%            'Columns'       - Columns to retrieve. Default is:
%                          {'ALPHA_J2000','DELTA_J2000','X_IMAGE','Y_IMAGE','MAG_APER1','MAG_APER2','MAGERR_APER1','MAGERR_APER2','MU_MAX','A_WORLD','B_WORLD','THETA_WORLD','ELONGATION','CLASS_STAR'};
%            'Table'         - Table to search. Default is 'HSCv1'.
% Output : - An AstCat object with element per search table.
%            Each element contains the list of sources found in search
%            table.
%          - Cell array of query result string.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Out,QS]=VO.MAST.wget_hsc_sources(1,1)
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;


DefV.HalfSize             = 10;
DefV.HalfSizeUnits        = 'arcsec';
DefV.SearchTable          = {'Catalog_WFPC2_SourceExtractor','Catalog_ACS_SourceExtractor','Catalog_WFC3_SourceExtractor'};
DefV.Columns              = {'ALPHA_J2000','DELTA_J2000','X_IMAGE','Y_IMAGE','MAG_APER1','MAG_APER2','MAGERR_APER1','MAGERR_APER2','MU_MAX','A_WORLD','B_WORLD','THETA_WORLD','ELONGATION','CLASS_STAR'};
DefV.Table                = 'HSCv1';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

RA  = celestial.coo.convertdms(RA,'gH','r');  % [rad]
Dec = celestial.coo.convertdms(Dec,'gD','R'); % [rad]

HalfSize = convert.angular(InPar.HalfSizeUnits,'rad',InPar.HalfSize);

[RA1,RA2,Dec1,Dec2] = celestial.coo.coo2box(RA,Dec,HalfSize,'deg'); % output in deg

Ntable = numel(InPar.SearchTable);
Out    = AstCat(Ntable,1);
for Itable=1:1:Ntable
    Query{1} = InPar.Columns;
    Query{2} = InPar.SearchTable(Itable);
    Query{3} = sprintf('ALPHA_J2000>%f and ALPHA_J2000<%f and DELTA_J2000> %f and DELTA_J2000<%f',RA1,RA2,Dec1,Dec2);
    [Out(Itable),ColCell,Status,QS{Itable},QueryStr] = VO.MAST.query_casjobs(Query,'Table',InPar.Table);
    Out(Itable).ColCell = ColCell;
end

