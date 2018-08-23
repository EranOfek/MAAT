function [RefMatch,RefCat,Cat,SortedCat]=cat_match_para(RefCat,RefColCoo,Cat,ColCoo,SearchRad,ExcludeSelf,Nproc);
%----------------------------------------------------------------------------
% cat_match_para function                                          Catalogue
% Description: Match 2 catalogs by object positions - a parallel version.
%              See cat_match.m
% Input  : - Reference catalog
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
%          - Number of processors. Default is 7.
% Output : - Structure containing the matches:
%            .IndexCat  - index of nearst match in the sorted search catalog.
%            .Dist      - distance to the nearest match [radians].
%            .PA        - PA of the nearest match [radians].
%            .Nmatch    - number of matches.
%          - Original reference catalog.
%          - Sorted search catalog.
%          - Search catalog sorted according to the Reference catalog.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                   September 2009
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------------
RAD = 180./pi;

Def.ExcludeSelf = 0;
Def.Nproc       = 7;
if (nargin==5),
   ExcludeSelf = Def.ExcludeSelf;
   Nproc       = Def.Nproc;
elseif (nargin==6),
   Nproc       = Def.Nproc;
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end

% sort by declination
%Cat     = sortrows(Cat,ColCoo(2));

NRefCat = size(RefCat,1);
NCat    = size(Cat,1);


RefMatchIndexCat = zeros(NRefCat,1).*NaN;
RefMatchDist     = zeros(NRefCat,1).*NaN;
RefMatchPA       = zeros(NRefCat,1).*NaN;
RefMatchNmatch   = zeros(NRefCat,1);
matlabpool(Nproc);
parfor I=1:1:NRefCat,
   [Lines,DistL,PAL]=cat_search(Cat,ColCoo,[RefCat(I,RefColCoo)],SearchRad);
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
   if (length(Lines)==0),
      % no match
      RefMatchNmatch(I) = 0;
   elseif (length(Lines)==1),
      % single match
      RefMatchNmatch(I) = 1;
      RefMatchIndexCat(I) = Lines;
      RefMatchDist(I)     = DistL;
      RefMatchPA(I)       = PAL;
   else
      % multi match
      RefMatchNmatch(I) = length(Lines);
      [Min,MinI] = min(DistL);
      RefMatchIndexCat(I) = Lines(MinI);
      RefMatchDist(I)     = DistL(MinI);
      RefMatchPA(I)       = PAL(MinI);
   end
end

RefMatch.IndexCat = RefMatchIndexCat;
RefMatch.Dist     = RefMatchDist;
RefMatch.PA       = RefMatchPA;
RefMatch.Nmatch   = RefMatchNmatch;

matlabpool close;

%RefMatch.IndexCat
Inan = find(isnan(RefMatch.IndexCat)==1);
Temp.IndexCat = RefMatch.IndexCat;
Temp.IndexCat(Inan) = 1;
SortedCat = Cat(Temp.IndexCat,:);
SortedCat(Inan,:) = NaN;
