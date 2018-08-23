function [SpDist,SpPA]=cat_acorr(VecRA,VecDec,SearchRad);
%-------------------------------------------------------------------------------
% cat_acorr function                                                  Catalogue
% Description: Auto correlate a catalog with itself. Given a catalog, look
%              for all neighbors found to the North of each source, within
%              a given distance.
% Input  : - Column vector of RA [rad], (sorted by declination).
%          - Column vector of Dec [rad].
%          - Search radius [rad].
% Output : - Sparse matrix containing the distances [radians] between
%            all selected neighboors which their distance is smaller
%            than the search radius.
%            Only the upper triangle of the matrix is filled.
%          - Sparse matrix containing position angle between neighboors [rad].
% Tested : Matlab 7.3
%     By : Eran O. Ofek                January 2008
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-------------------------------------------------------------------------------



N = length(VecDec);

SpDist = sparse(N,N,ceil(0.1.*N));
SpPA   = sparse(N,N,ceil(0.1.*N));

for I=1:1:N-1,
   [L,Dist,PA] = cat_search([VecRA, VecDec],[1 2],[VecRA(I),VecDec(I)],SearchRad);
   If   = find(L>I);
   if (isempty(If)==0),
      Dist = Dist(If);
      PA   = PA(If);
      Js   = L(If);
      SpDist(I,Js) = Dist;
      SpPA(I,Js)   = PA;
   end
end




