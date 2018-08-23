function [OutRef,Extra]=xcat(Sim,varargin)
% Cross match an AstCat object with an external catalog.
% Package: @AstCat
% Description: Cross match an AstCat object with an external catalog.
%              For each source in an AstCat (or SIM) object, search for
%              source in an external catalog (e.g., SDSS).
%              Generate a an AstCat object in which for each source in the
%              input AstCat object the nearest source within the search
%              radius, in the external catalog is listed.
%              If no match, then line is populated with NaNs.
% Input  : - An AstCat (or SIM) object. Each element in the object contains
%            a catalog, and each catalog will be cross matched.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'RA' - J2000.0 R.A. [radians] around which to query the
%                   external catalog. If empty, then will attempt to obtain
%                   the center of the catalog from image header.
%                   Default is empty.
%            'Dec' - J2000.0 Dec. [radians] around which to query the
%                   external catalog. If empty, then will attempt to obtain
%                   the center of the catalog from image header.
%                   Default is empty.
%            'Radius' - Extraction radius for external catalog [radians].
%                   If empty, then will attempt to obtain
%                   the center of the catalog from image header.
%                   Default is empty.
%            'CCDSEC' - CCDSEC header keyword. If empty, then use actual
%                   image size. Default is 'CCDSEC'.
%            'Refcat' - External catalog: 'sdss' | 'apass' | ...
%                   See VO.search.get_cat for options.
%                   Default is 'sdss'.
%            'SearchRad' - Match radius. Default is 2.
%            'SearchRadUnits' - Match radius units. Default is 'arcsec'.
%            'GenCatFun' - If catalog is empty then try to generate catalog
%                   using this function. Default is @mextractor.
%            'GenCatPar' - Cell array of parameters to pass to the catalog
%                   extraction program. Default is {}.
%            'CatColRA' - Cell array of column names for the J2000.0 R.A.
%                   in the input catalog. The program will search for the
%                   first existing column name.
%                   Default is {'ALPHAWIN_J2000','RA'}.
%            'CatColDec' - Like 'CatColRA', but for the J2000.0 Dec.
%                   Default is {'DELTAWIN_J2000','Dec'}.
%            'CatUnits' - Units of the coordinates in the input catalog.
%                   Default is 'deg'.
%            'RefColRA' - Cell array of column names for the J2000.0 R.A.
%                   in the external catalog. The program will search for 
%                   the first existing column name.
%                   Default is {'RA'}.
%            'RefColDec' - Like 'RefColRA', but for the J2000.0 Dec.
%                   Default is {'Dec'}.
%            'RefUnits' - Units of the coordinates in the external catalog.
%                   Default is 'rad'.
%            'EmptyVal' - Value of non matched sources. Default is NaN.
% Output : - An AstCat object. Each element corresponds to an element in
%            the input AstCat object. Each element contains a acatalog in
%            which for each source in the input catalog the external
%            catalog entry is provided. NaNs if no match.
%          - A structure array of extra information.
%            Each element corresponds to an element in the input AstCat
%            object. The following fields are available:
%            'Nmatch' - number of matches for each source.
%            'Dist' - angular distance to nearest match [radians].
%            'PA' - position angle for nearest match [radians].
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [OutRef,Extra]=xcat(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

%AstC = sim_xcat(AstC,varargin{:});

CatField     = AstCat.CatField;
ColField     = AstCat.ColField;
ColCellField = AstCat.ColCellField;


DefV.RA                   = [];
DefV.Dec                  = [];
DefV.Radius               = [];
DefV.KeyCCDSEC            = 'CCDSEC';   % if empty use actual image size
DefV.RefCat               = 'SDSSDR10';     % AstCat or catalaog name
DefV.SearchRad            = 2;
DefV.SearchRadUnits       = 'arcsec';
% generate catalog
DefV.GenCatFun            = @mextractor;
DefV.GenCatPar            = {};
% catalaog attributes
DefV.CatColRA             = {'ALPHAWIN_J2000','RA'};
DefV.CatColDec            = {'DELTAWIN_J2000','Dec'};
DefV.CatUnits             = 'deg';
% xcat
DefV.RefColRA             = {'RA'};
DefV.RefColDec            = {'Dec'};
DefV.RefUnits             = 'rad';
DefV.EmptyVal             = NaN;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

InPar.SearchRad = convert.angular(InPar.SearchRadUnits,'rad',InPar.SearchRad);

ConvertCat = convert.angular(InPar.CatUnits,'rad');
ConvertRef = convert.angular(InPar.RefUnits,'rad');

OutRef = AstCat(size(Sim));
Extra  = Util.struct.struct_def({'Nmatch','Dist','PA'},size(Sim));
Nsim = numel(Sim);
for Isim=1:1:Nsim
    % get RA/Dec/Radius from image information
    if (isempty(InPar.RA) || isempty(InPar.Dec) || isempty(InPar.Radius))
        Foot = footprint(Sim(Isim),'CCDSEC',InPar.KeyCCDSEC,'OutUnits','rad');
        InPar.RA      = Foot.GcenLong;
        InPar.Dec     = Foot.GcenLat;
        InPar.Radius  = Foot.Radius;
    end
    
    % check if catalog exist - if not generate catalog
    if isempty(Sim(Isim).(CatField))
        % generate catalog
        Sim(Isim) = InPar.GenCatFun(Sim(Isim),InPar.GenCatPar{:});
    end
    
    
    
    % get relevant columns from catalog
    [~,CatColRA]     = select_exist_colnames(Sim,InPar.CatColRA');
    [~,CatColDec]    = select_exist_colnames(Sim,InPar.CatColDec');
    

    % get reference catalog
    if (ischar(InPar.RefCat))
        %Ref = VO.search.get_cat(InPar.RefCat,InPar.RA,InPar.Dec,InPar.Radius,'RadUnits','rad');
        Ref = VO.search.cat_cone(InPar.RefCat,InPar.RA,InPar.Dec,InPar.Radius,'RadiusUnits','rad','OutType','astcat');
    else
        error('RefCat type not supported');
    end
    
    % get relevant columns from Reference
    [~,RefColRA]     = select_exist_colnames(Ref,InPar.RefColRA');
    [~,RefColDec]    = select_exist_colnames(Ref,InPar.RefColDec');
    
    % sort catalogs
    if (isempty(Ref.(CatField)))
        error('Ref catalog is empty');
    end
    Ref = sortrows(Ref,RefColDec);
    
    
    [Ind,Extra(Isim).Nmatch,Extra(Isim).Dist,Extra(Isim).PA] = celestial.search.match_coo_nearest(Sim(Isim).(CatField)(:,[CatColRA, CatColDec]).*ConvertCat,...
                                                           Ref.(CatField)(:,[RefColRA, RefColDec]).*ConvertRef,...
                                                           InPar.SearchRad);
    
    
    Nsrc = size(Sim(Isim).(CatField),1);
    Ncol = size(Ref.(CatField),2);
    SrcInd = (1:1:Nsrc).';
    NN = ~isnan(Ind);
    
    
    if isnan(InPar.EmptyVal)
        OutRef(Isim).(CatField) = nan(Nsrc,Ncol);
    else
        OutRef(Isim).(CatField) = InPar.EmptyVal.*ones(Nsrc,Ncol);
    end
    
    OutRef(Isim).(CatField)(SrcInd(NN),:) = Ref.(CatField)(Ind(NN),:);
    OutRef(Isim).(ColField) = Ref.(ColField);
    OutRef(Isim).(ColCellField) = Ref.(ColCellField);
end
                    
                                     