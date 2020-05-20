function [Ind,FlagUnique]=search_sortedY_multi(Cat,Long,Lat,Radius,FlagUnique)
% Search a single X/Y in a catalog sorted by Y (planar geometry)
% Package: VO.search
% Description: A low level function for a single cone search
%              in a [X, Y] array.
% Input  : - An array of [X, Y] in radians, sorted by Y.
%            The program doesnot verify that the array is sorted.
%          - X to search.
%          - Y to search.
%          - Search radius.
%            If radius is negative then will add .Dist to the output index
%            structure.
%          - A vector of false with the length of the catalog,
%            after first iteration provide the output.
% Output : - A strucure array in which the number of elements equal to the
%            number of searched coordinates, and with the following
%            fields:
%            'Ind' - Vector of indices of matched sources.
%            'Nmatch' - Number of matched sources.
%            'Dist' - Distance between sources. This is provided only if
%                   the input search radius is negative.
%          - A logical vector of length equal to the number of searched
%            coordinates (Y) which flag the first unique source.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=sortrows(rand(10000,2),2);
%          Ind=VO.search.search_sortedY_multi(Cat,0.5,0.5,0.01)
%          Ind=VO.search.search_sortedY_multi(Cat,Cat(:,1),Cat(:,2),0.0001)
% Reliable: 2


if Radius<0
    CalcDist = true;
    Radius = abs(Radius);
else
    CalcDist = false;
end

if (nargin<5)
    FlagUnique = false(size(Cat,1),1);
end

Radius2 = Radius.^2;

Col.Lon = 1;
Col.Lat = 2;

% somewhat slower version:
% Ilow  = Util.find.bin_sear(Cat(:,Col.Lat),Lat-Radius);
% Ihigh = Util.find.bin_sear(Cat(:,Col.Lat),Lat+Radius);

Lat   = Lat(:).';  % convert Lat to a row vector
Nlat  = numel(Lat); % number of latitudes to search
Ilat  = [(1:1:Nlat).', (1:1:Nlat).'+Nlat];

Ncat  = size(Cat,1);
Inear = Util.find.mfind_bin(Cat(:,Col.Lat),[Lat-Radius, Lat+Radius]);

% Inear(Ilat) is a two column matrix [low, high] index for each latitud
% search
Ilowhigh = double(Inear(Ilat));
Ilow     = Ilowhigh(:,1);
Ihigh    = min(Ncat,Ilowhigh(:,2)+1); % add 1 because of the way mfind_bin works

Ind = Util.struct.struct_def({'Ind','Nmatch'},Nlat,1);
for I=1:1:Nlat
    %Dist = celestial.coo.sphere_dist_fast(Long(I),Lat(I), Cat(Ilow(I):Ihigh(I),Col.Lon), Cat(Ilow(I):Ihigh(I),Col.Lat));
    Dist2 = (Long(I) - Cat(Ilow(I):Ihigh(I),Col.Lon)).^2 + (Lat(I) - Cat(Ilow(I):Ihigh(I),Col.Lat)).^2;
    FlagDist = Dist2 <= Radius2;
    Ind(I).Ind    = Ilow(I)-1+find(FlagDist);
    %Ind(I).Ind    = Ilow(I)-1+find(Dist <= Radius);
    Ind(I).Nmatch = numel(Ind(I).Ind);
    if CalcDist
        Ind(I).Dist   = sqrt(Dist2(FlagDist));
    end
    
    if ~any(FlagUnique(Ind(I).Ind))
        % a new unique source
        FlagUnique(I) = true;
    end
    
end