function [IndTable,CatFlagNearest,CatFlagAll]=search_sortedlat_multiNearest(Cat,Long,Lat,Radius,DistFun)
% Search a single long/lat in a catalog sorted by latitude
% Package: VO.search
% Description: A low level function for a single cone search
%              in a [Long, Lat] array.
% Input  : - An array of [Long, Lat] in radians, sorted by Lat.
%            The program doesnot verify that the array is sorted.
%          - Longitude [radians] to search.
%          - Latitude [radians] to search.
%          - Radius [radians] to search.
%          - A function handle for calculating distances Fun(X1,Y1,X2,Y2).
%            Default is @celestial.coo.sphere_dist_fast.
% Output : - A three column matrix with, one line per line in Long,Lat.
%            Columns are [Index of nearest source, within search radius, in Cat,
%            Distance, Total number of matches within radius].
%          - A vector of logical (length as Cat), which indicate the object
%            in Cat that were identified as the nearest object to a source in
%            Long, Lat.
%          - The same as the previous output, but for all the sources
%            within the search radius.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=sortrows(rand(10000,2),2);
%          Ind=VO.search.search_sortedlat_multi(Cat,0.5,0.5,0.01)
% Reliable: 2

arguments
    Cat
    Long
    Lat
    Radius
    DistFun function_handle      = @celestial.coo.sphere_dist_fast;
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

IndTable = [nan(Nlat,2), zeros(Nlat,1)]; % [Index, Dist, Nmatch]
%Util.struct.struct_def({'Ind','Nmatch','Dist'},Nlat,1);

CatFlagNearest  = false(Ncat,1);
CatFlagAll      = false(Ncat,1);


for I=1:1:Nlat
    %Dist  = celestial.coo.sphere_dist_fast(Long(I),Lat(I), Cat(Ilow(I):Ihigh(I),Col.Lon), Cat(Ilow(I):Ihigh(I),Col.Lat));
    Dist  = DistFun(Long(I),Lat(I), Cat(Ilow(I):Ihigh(I),Col.Lon), Cat(Ilow(I):Ihigh(I),Col.Lat));
    FlagDist = Dist <= Radius;
    
    IndI  = Ilow(I)-1+find(FlagDist);
    DistI = Dist(FlagDist);
    if ~isempty(DistI)
        [MinDist, MinInd] = min(DistI);
        IndTable(I,1) = IndI(MinInd);   % VERIFY THIS???
        IndTable(I,2) = MinDist;
        IndTable(I,3) = numel(IndI);
        
        CatFlagNearest(IndTable(I,1)) = true;
        CatFlagAll(IndI)              = true;
        
    end
    
end
