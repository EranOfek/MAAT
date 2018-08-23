function Ind=search_sortedlat(Cat,Long,Lat,Radius)
% Search a single long/lat in a catalog sorted by latitude
% Package: VO.search
% Description: A low level function for a single cone search
%              in a [Long, Lat] array.
% Input  : - An array of [Long, Lat] in radians, sorted by Lat.
%            The program doesnot verify that the array is sorted.
%          - Longitude [radians] to search.
%          - Latitude [radians] to search.
%          - Radius [radians] to search.
% Output : - Indices of the entries in the input [Lon, Lat] catalog which
%            are found within Radius of Long,Lat.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=sortrows(rand(10000,2),2);
%          Ind=VO.search.search_sortedlat(Cat,0.5,0.5,0.01)
% Reliable: 2

Col.Lon = 1;
Col.Lat = 2;

% somewhat slower version:
% Ilow  = Util.find.bin_sear(Cat(:,Col.Lat),Lat-Radius);
% Ihigh = Util.find.bin_sear(Cat(:,Col.Lat),Lat+Radius);

Ncat  = size(Cat,1);
Inear = Util.find.mfind_bin(Cat(:,Col.Lat),[Lat-Radius, Lat+Radius]);
Ilow  = double(Inear(1));
Ihigh = min(Ncat,double(Inear(2)) + 1);  % add 1 because of the way mfind_bin works

Dist  = celestial.coo.sphere_dist_fast(Long,Lat, Cat(Ilow:Ihigh,Col.Lon), Cat(Ilow:Ihigh,Col.Lat));

Ind   = Ilow-1+find(Dist <= Radius);





