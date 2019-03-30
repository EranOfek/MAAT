function [AstC,ColCell]=cat_cone(CatName,RA,Dec,Radius,varargin)
% Cone search a local or online catalog.
% Package: VO.search
% Description: A uniform interface function for catalogs cone search.
%              The catalog can be either a local catalog in the +cats
%              package, or an online (web) catalog for which inteface
%              function exist, or a local catalog in HDF5/HTM format.
% Input  : - catalog name, or function handle.
%            The following options are available:
%            (1) To access catalogs in the +cats directory use e.g.,
%                @cats.X.ROSAT_faint.
%            (2) To access catalogs which have access program use
%                function handle.
%            (3) To access catalogs in catsHTM format, specify catalog
%                name.
%          - J2000.0 R.A. [radians, [H M S], or sexagesimal string].
%          - J2000.0 Dec. [radians, [sign D M S], or sexagesimal string].
%          - Search radius [arcsec].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'RadiusUnits'  - Search radius units. Default is 'arcsec'.
%            'OutType'      - 'mat' | 'astcat'. Default is 'mat'.
% Output : - The catalog.
%          - Cell array of column names.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cat,ColCell] = VO.search.cat_cone(@cats.X.ROSAT_faint,1,1,1)
%          [Cat,ColCell] = VO.search.cat_cone(@cats.X.ROSAT_faint,1,1,1,'OutType','astcat')
%          [Cat,ColCell] = VO.search.cat_cone('UCAC4',1,1,1)
%          [Cat,ColCell] = VO.search.cat_cone('UCAC4',1,1,1,'OutType','astcat')
%          [Cat,ColCell] = VO.search.cat_cone('GAIADR1',1,1,1,'OutType','astcat')
% Reliable: 2
%--------------------------------------------------------------------------

CatField     = AstCat.CatField;
ColCellField = AstCat.ColCellField;
ColUnitsField= AstCat.ColUnitsField;


DefV.RadiusUnits          = 'arcsec';
DefV.OutType              = 'mat';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%FFU
% deal with box searches
%BoxSize   = Radius;
%Radius    = sqrt(sum(BoxSize.^2));

RadiusAS  = convert.angular(InPar.RadiusUnits,'arcsec',Radius);

if (isa(CatName,'function_handle'))
    % function handle
    AstC = CatName(RA,Dec,RadiusAS);
    % default output is AstCat
    ColCell = AstC.(ColCellField);
    switch lower(InPar.OutType)
        case 'mat'
            ColCell = AstC.(ColCellField);
            AstC    = AstC.(CatField);
        otherwise
            % do nothing
    end
else
    % if catalog is a string than assume this is an HDF5/HTM catalog
    [Cat,ColCell,ColUnits] = catsHTM.cone_search(CatName,RA,Dec,RadiusAS);
    % default output is matrix
    switch lower(InPar.OutType)
        case 'astcat'
            AstC = AstCat;
            AstC.(CatField)     = Cat;
            AstC.(ColCellField) = ColCell;
            AstC.(ColUnitsField)= ColUnits;
            AstC = colcell2col(AstC);
        otherwise
            AstC = Cat; 
    end
end

