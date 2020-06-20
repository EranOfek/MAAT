function [Links,GALEX]=searchGALEXimg(RA,Dec,varargin)
% Search GALEX images by coordinates
% Package: VO.GALEX
% Description: Search GALEX images by coordinates in the
%              cats.GALEX.GALEXimg local catalog.
% Input  : - J2000.0 R.A. [rad, [H M S], or sexagesimal
%            string], or and Index in the GALEX images file.
%          - J2000.0 Dec. [rad, [Sign D M S], or sexagesimal
%            string]. If first argument is ID then this need to
%            be an empty matrix. Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SearchRadius' - Default is 1.
%            'SearchRadiusUnits' - Default is 'deg';
%            'CooUnits' - {'deg' | 'rad'}, Default is 'rad'.
%                         If 'rad' then will use celestial.coo.convertdms
%                         to read coordinates.
% Output : - A structure of links to the various products.
%          - An AstCat object containing the selected images.
% Example: [Links,GALEX]=VO.GALEX.searchGALEXimg(RA,Dec)
% Reliable: 2


RAD = 180./pi;
BaseURL = 'http://galex.stsci.edu/data/';


if (nargin==1)
    Dec = [];
end

DefV.SearchRadius       = 1;
DefV.SearchRadiusUnits  = 'deg';
DefV.CooUnits           = 'rad';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


switch lower(InPar.CooUnits)
    case 'deg'
        RA  = RA./RAD;
        Dec = Dec./RAD;
    otherwise
        RA  = celestial.coo.convertdms(RA,'gH','r');
        Dec = celestial.coo.convertdms(Dec,'gD','R');
end

InPar.SearchRadius = convert.angular(InPar.SearchRadiusUnits,'rad',InPar.SearchRadius);
GALEXimg = cats.GALEX.GALEXimg;
Ind = VO.search.search_sortedlat([GALEXimg.Cat.RA, GALEXimg.Cat.Dec],RA,Dec,InPar.SearchRadius);

GALEX = GALEXimg;
GALEX.Cat = GALEX.Cat(Ind,:);
Nind = numel(Ind);
Links.intNUV = cell(Nind,1);
Links.intFUV = cell(Nind,1);
for Iind=1:1:Nind
    Links.intNUV{Iind} = sprintf('%s%s',BaseURL,GALEX.Cat.nuv_fileNPath{Iind});
    Links.intFUV{Iind} = sprintf('%s%s',BaseURL,GALEX.Cat.fuv_fileNPath{Iind});
end

% create links to other file types:
Links.cntNUV  = regexprep(Links.intNUV,'-int.fits','-cnt.fits');
Links.cntFUV  = regexprep(Links.intFUV,'-int.fits','-cnt.fits');
Links.rrhrNUV = regexprep(Links.intNUV,'-int.fits','-rrhr.fits');
Links.rrhrFUV = regexprep(Links.intFUV,'-int.fits','-rrhr.fits');
Links.mcat    = regexprep(Links.intNUV,'nd-int.fits','-xd-mcat.fits');

% clean the nulls
NullStr = sprintf('%s%s',BaseURL,'null');
Links.intNUV  = regexprep(Links.intNUV,NullStr,'');
Links.intFUV  = regexprep(Links.intFUV,NullStr,'');
Links.cntNUV  = regexprep(Links.cntNUV,NullStr,'');
Links.cntFUV  = regexprep(Links.cntFUV,NullStr,'');
Links.rrhrNUV = regexprep(Links.rrhrNUV,NullStr,'');
Links.rrhrFUV = regexprep(Links.rrhrFUV,NullStr,'');
Links.mcat    = regexprep(Links.mcat,NullStr,'');




