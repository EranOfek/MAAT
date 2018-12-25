function AM=match_with_cat(AC,CatName,varargin)
% Match an AstCat object against an catsHTM catalog
% Package: @AstCat
% Description: Match a source catalog (with RA/Dec) stored in an AstCat
%              object with an external catalog in the catsHTM format.
%              For each input AstCat element, the code return an AstCat
%              element with the same number of lines (i.e., sources) with
%              the information from the catsHTM catalog, and some metadata
%              (e.g., distance, number of matches).
% Input  : - An AstCat object. Each element contains a catalog with RA/Dec
%            columns.
%          - catsHTM catalog name (e.g., 'GAIADR2').
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SearchRadius' - Match search radius. Default is 3.
%            'SearchRadiusUnits' - Match search radius units. Default is
%                            'arcsec'.
%            'ColRA' - Cell array of possible J2000.0 R.A. column names in
%                      the input AstCat object.
%            'ColDec' - Cell array of possible J2000.0 Dec. column names in
%                      the input AstCat object.
%            'CooUnits' - Coordinates units in the input AstCat object.
%                      Default is 'deg'.
% Output : - An AstCat object with the same number of elements as the input
%            AstCat object. Each element contains the same number of lines
%            (i.e., sources) as the input object, but with data from the
%            catsHTM catalog.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AM=match_with_cat(S(1),'GAIADR2');
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;
ARCSEC_DEG = 3600;

CatField     = 'Cat';
ColField     = 'Col';
ColCellField = 'ColCell';

DefV.SearchRadius         = 3;
DefV.SearchRadiusUnits    = 'arcsec';
DefV.ColRA                = {'ALPHAWIN_J2000','RA','ALPHA_J2000'};
DefV.ColDec               = {'DELTAWIN_J2000','Dec','DELTA_J2000'};
DefV.CooUnits             = 'deg';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


InPar.SearchRadius = convert.angular(InPar.SearchRadiusUnits,'arcsec',InPar.SearchRadius,InPar.SearchRadiusUnits);

Nac = numel(AC);
% for each AstCat element
for Iac=1:1:Nac
    
    [~,ColRA,~]  = select_exist_colnames(AC(Iac),InPar.ColRA');
    [~,ColDec,~] = select_exist_colnames(AC(Iac),InPar.ColDec');
    
    RA  = convert.angular(InPar.CooUnits,'rad',AC(Iac).(CatField)(:,ColRA));
    Dec = convert.angular(InPar.CooUnits,'rad',AC(Iac).(CatField)(:,ColDec));
    
    % define search region
    % RA/Dec to cosine directions
    [CD1,CD2,CD3] = celestial.coo.coo2cosined(RA, Dec);
    % approximate catalog center
    [MeanRA,MeanDec] = celestial.coo.cosined2coo(nanmedian(CD1),nanmedian(CD2),nanmedian(CD3));
    Dist             = celestial.coo.sphere_dist(MeanRA,MeanDec,...
                                                 RA, Dec);
    % search radius for possible matching sources
    BigRadius = max(Dist);
    
    % get entire catalog overlap to the AstCat object
    % ccordinates units are 'rad'
    [HCat,Col] = catsHTM.cone_search(CatName,MeanRA,MeanDec,BigRadius.*RAD.*ARCSEC_DEG,'OutType','astcat');
    HCat       = sortrows(HCat,2);
    
    CatColNames = [InPar.ColRA(:), InPar.ColDec(:)];
    [CatUnits{1:numel(InPar.ColRA)}] = deal(InPar.CooUnits);
    [CooType{1:numel(InPar.ColRA)}] = deal('sphere');
    
    % find a match to each HCat source - output is the size of HCat
    %[AM] = match(AC(Iac),HCat,'SkipWCS',false,'CatColNames',CatColNames,'CatUnits',CatUnits,'CooType',CooType);
    % find a match to each AC source - output is the size of AC(Iac)
    
    [AM(Iac)] = match(HCat,AC(Iac),'SkipWCS',false,'CatColNames',CatColNames,...
                                   'RefUnits',CatUnits,'CatUnits',{'rad','rad','rad'}','CooType',CooType',...
                                   'SearchRad',InPar.SearchRadius,'SearchRadUnits',InPar.SearchRadiusUnits);
      
    
end
    
    