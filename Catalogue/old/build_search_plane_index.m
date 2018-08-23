function [CellCat,IndCat]=build_search_plane_index(Catalog,ColCoo,VecXb,VecYb);

CellCat.Catalog = Catalog;
Nc              = size(Catalog,1);
CellCat.Ind     = zeros(Nc,1);
Nxb = length(VecXb);
Nyb = length(VecYb);
IndCat = zeros((Nxb-1).*(Nyb-1),2);
K = 0;
for Ixb=1:1:Nxb-1,
   for Iyb=1:1:Nyb-1,
      K = K + 1;
      Ic = find(Catalog(:,ColCoo(1))>=VecXb(Ixb) & ...
                Catalog(:,ColCoo(1))<VecXb(Ixb+1) & ...
                Catalog(:,ColCoo(2))>=VecYb(Iyb) & ...
                Catalog(:,ColCoo(2))<VecYb(Iyb+1));

      CellCat.Ind(Ic) = K;
      IndCat(K,1) = (VecXb(Ixb) + VecXb(Ixb+1)).*0.5;
      IndCat(K,2) = (VecYb(Iyb) + VecYb(Iyb+1)).*0.5;
   end
end





