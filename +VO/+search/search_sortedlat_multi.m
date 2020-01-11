function [Ind,FlagUnique]=search_sortedlat_multi(Cat,Long,Lat,Radius,FlagUnique)
% Search a single long/lat in a catalog sorted by latitude
% Package: VO.search
% Description: A low level function for a single cone search
%              in a [Long, Lat] array.
% Input  : - An array of [Long, Lat] in radians, sorted by Lat.
%            The program doesnot verify that the array is sorted.
%          - Longitude [radians] to search.
%          - Latitude [radians] to search.
%          - Radius [radians] to search.
%          - A vector of false with the length of the catalog,
%            after first iteration provide the output.
% Output : - A strucure array in which the number of elements equal to the
%            number of searched coordinates (Lat), and with the following
%            fields:
%            'Ind' - Vector of indices of matched sources.
%            'Nmatch' - Number of matched sources.
%          - A logical vector of length equal to the number of searched
%            coordinates (Lat) which flag the first unique source.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=sortrows(rand(10000,2),2);
%          Ind=VO.search.search_sortedlat(Cat,0.5,0.5,0.01)
% Reliable: 2

if (nargin<5)
    FlagUnique = false(size(Cat,1),1);
end

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
    Dist  = celestial.coo.sphere_dist_fast(Long(I),Lat(I), Cat(Ilow(I):Ihigh(I),Col.Lon), Cat(Ilow(I):Ihigh(I),Col.Lat));
    Ind(I).Ind    = Ilow(I)-1+find(Dist <= Radius);
    Ind(I).Nmatch = numel(Ind(I).Ind);
    if ~any(FlagUnique(Ind(I).Ind))
        % a new unique source
        FlagUnique(I) = true;
    end
    
end