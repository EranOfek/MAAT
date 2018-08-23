function [RefMatch,RefCat,Cat,SortedCat,CatNoMatchInd]=cat_match(RefCat,RefColCoo,Cat,ColCoo,SearchRad,ExcludeSelf,SearchShape,WorkCoo,CooType)
%------------------------------------------------------------------------------
% cat_match function                                                 Catalogue
% Description: Match two catalogs by object positions. If the catalogs are
%              not sorted then the program will sort them.
% Input  : - Reference catalog.
%          - Columns of [RA, Dec] in the reference catalog [rad].
%          - Search catalog sorted by declination.
%          - Columns of [RA, Dec] in the search catalog [rad].
%          - Search radius [radians]
%          - Exclude self sources {2 | 1 | 0}.
%            0 - do not exclude source (default).
%            1 - exclude sources found within 1e-10 rad
%                 from the searched source (This is for comparison
%   	          of a catalog with itself).
%            2 - exclude the source itself (by its index).
%          - Shape of search region: {'circle' | 'box'},
%            default is 'circle'.
%          - Coordinate to search by {'RA' | 'Dec'}, default is 'Dec'.
%          - Coordinates type:
%            'sphere' - Spherical coordinates. Coordinates must be
%                       given in radians. Default.
%            'plane'  - plan coordinates.
% Output : - Structure containing the matches as found in the sorted catalogs.
%            The structure contains the following fields:
%            .IndexCat  - For each source in the sorted version of RefCat
%                         this is the index of the nearest source in Cat
%                         found within the search radius.
%                         Therefore RefCat(1,:) is matched to
%                         RefCat(RefMatch.IndexCat(1)).
%                         If the source was not found in Cat, then the
%                         index is NaN.
%            .Dist      - Distance to the nearest match [radians].
%            .PA        - PA of the nearest match [radians].
%            .Nmatch    - Number of matches.
%          - Sorted reference catalog.
%          - Sorted search catalog.
%          - Search catalog sorted according to the Reference catalog.
%          - Indices of sources in Cat which does not have a match in
%            RefCat.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Sep 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;

Def.ExcludeSelf = 0;
Def.SearchShape = 'circle';
Def.WorkCoo     = 'Dec';
Def.CooType     = 'sphere';

if (nargin==5),
   ExcludeSelf = Def.ExcludeSelf;
   SearchShape = Def.SearchShape;
   WorkCoo    = Def.WorkCoo;
   CooType     = Def.CooType;
elseif (nargin==6),
   SearchShape = Def.SearchShape;
   WorkCoo    = Def.WorkCoo;
   CooType     = Def.CooType;
elseif (nargin==7),
   WorkCoo    = Def.WorkCoo;
   CooType     = Def.CooType;
elseif (nargin==8),
   CooType     = Def.CooType;
elseif (nargin==9),
   % do nothing
else
   error('Illegal number of input arguments');
end

% sort by declination
%Cat     = sortrows(Cat,ColCoo(2));
switch lower(WorkCoo)
 case 'dec'
    IndWorkCoo = 2;
 case 'ra'
    IndWorkCoo = 1;
 otherwise
    error('Unknwon WorkCoo option');
end

if (~issorted(Cat(:,ColCoo(IndWorkCoo)))),
    % sort input catalog
    %error('Input Cat is not sorted');
    Cat     = sortrows(Cat,ColCoo(IndWorkCoo));
end

if (~issorted(RefCat(:,RefColCoo(IndWorkCoo)))),
    % sort input catalog
    %error('Input RefCat is not sorted');
    RefCat     = sortrows(RefCat,RefColCoo(IndWorkCoo));
end

NRefCat = size(RefCat,1);
NCat    = size(Cat,1);


RefMatch.IndexCat = zeros(NRefCat,1).*NaN;
RefMatch.Dist     = zeros(NRefCat,1).*NaN;
RefMatch.PA       = zeros(NRefCat,1).*NaN;
RefMatch.Nmatch   = zeros(NRefCat,1);
for I=1:1:NRefCat,
   [Lines,DistL,PAL]=cat_search(Cat,ColCoo,RefCat(I,RefColCoo),SearchRad,SearchShape,WorkCoo,CooType);
   switch ExcludeSelf
    case 2
       % exclude by object index
       Lines = Lines(Lines~=I);
       DistL = DistL(Lines~=I);
       PAL   = PAL(Lines~=I);
    case 1
       % exclude by coordinates
       Iok = find(DistL>1e-10);
       Lines = Lines(Iok);
       DistL = DistL(Iok);
       PAL   = PAL(Iok);
    otherwise
       % do nothing
   end
   if (isempty(Lines)),
      % no match
      RefMatch.Nmatch(I) = 0;
   elseif (length(Lines)==1),
      % single match
      RefMatch.Nmatch(I) = 1;
      RefMatch.IndexCat(I) = Lines;
      RefMatch.Dist(I)     = DistL;
      RefMatch.PA(I)       = PAL;
   else
      % multi match
      RefMatch.Nmatch(I) = length(Lines);
      [~,MinI] = min(DistL);
      RefMatch.IndexCat(I) = Lines(MinI);
      RefMatch.Dist(I)     = DistL(MinI);
      RefMatch.PA(I)       = PAL(MinI);
   end
end

%RefMatch.IndexCat
Inan = find(isnan(RefMatch.IndexCat));
Temp.IndexCat = RefMatch.IndexCat;
Temp.IndexCat(Inan) = 1;
SortedCat = Cat(Temp.IndexCat,:);
SortedCat(Inan,:) = NaN;

if (nargout>4),
   CatNoMatchInd = setdiff((1:1:NCat).',RefMatch.IndexCat);
end
