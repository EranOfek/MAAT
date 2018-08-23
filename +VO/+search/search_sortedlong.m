function Ind=search_sortedlong(Cat,Long,Lat,Radius)
% Search a single long/lat in a catalog sorted by longitude
% Package: VO.search
% Description: A low level function for a single cone search
%              in a [Long, Lat] array.
% Input  : - An array of [Long, Lat] in radians, sorted by Long.
%            The program doesnot verify that the array is sorted.
%          - Longitude [radians] to search.
%          - Latitude [radians] to search.
%          - Radius [radians] to search.
% Output : - Indices of the entries in the input [Lon, Lat] catalog which
%            are found within Radius of Long,Lat.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=sortrows(rand(10000,2),1);
%          Ind=VO.search.search_sortedlong(Cat,0.5,0.5,0.01)
% Reliable: 2

Col.Lon = 1;
Col.Lat = 2;

% somewhat slower version:
% Ilow  = Util.find.bin_sear(Cat(:,Col.Lat),Lat-Radius);
% Ihigh = Util.find.bin_sear(Cat(:,Col.Lat),Lat+Radius);

Ncat  = size(Cat,1);

MaxLat = min(abs(Lat)+Radius,pi./2);
UpdatedRadius = Radius./cos(MaxLat);

if (UpdatedRadius>pi)
    % no need to search - use all sources
    VecI  = (1:1:Ncat).';
    Dist  = celestial.coo.sphere_dist_fast(Long,Lat, Cat(VecI,Col.Lon), Cat(VecI,Col.Lat));
    Ind   = find(Dist <= Radius);
else
    MaxLong = Long + UpdatedRadius;
    MinLong = Long - UpdatedRadius;
    if (MaxLong>(2.*pi))
        MinLong1 = 0;
        MaxLong1 = mod(MaxLong,2.*pi);
        
        MinLong2 = MinLong;
        MaxLong2 = 2.*pi;
        Inear1 = Util.find.mfind_bin(Cat(:,Col.Lon),[MinLong1, MaxLong1]);
        Inear2 = Util.find.mfind_bin(Cat(:,Col.Lon),[MinLong2, MaxLong2]);
        Ilow1  = double(Inear1(1));
        Ihigh1 = min(Ncat,double(Inear1(2)) + 1);  % add 1 because of the way mfind_bin works
        Ilow2  = double(Inear2(1));
        Ihigh2 = min(Ncat,double(Inear2(2)) + 1);  % add 1 because of the way mfind_bin works
        Dist1  = celestial.coo.sphere_dist_fast(Long,Lat, Cat(Ilow1:Ihigh1,Col.Lon), Cat(Ilow1:Ihigh1,Col.Lat));
        Dist2  = celestial.coo.sphere_dist_fast(Long,Lat, Cat(Ilow2:Ihigh2,Col.Lon), Cat(Ilow2:Ihigh2,Col.Lat));
        Ind    = [Ilow1-1+find(Dist1 <= Radius); Ilow2-1+find(Dist2 <= Radius)];
        
    elseif (MinLong<0)
        MinLong1 = 0;
        MaxLong1 = MaxLong;
        
        MinLong2 = mod(MinLong,2.*pi);
        MaxLong2 = 2.*pi;
        Inear1 = Util.find.mfind_bin(Cat(:,Col.Lon),[MinLong1, MaxLong1]);
        Inear2 = Util.find.mfind_bin(Cat(:,Col.Lon),[MinLong2, MaxLong2]);
        Ilow1  = double(Inear1(1));
        Ihigh1 = min(Ncat,double(Inear1(2)) + 1);  % add 1 because of the way mfind_bin works
        Ilow2  = double(Inear2(1));
        Ihigh2 = min(Ncat,double(Inear2(2)) + 1);  % add 1 because of the way mfind_bin works
        Dist1  = celestial.coo.sphere_dist_fast(Long,Lat, Cat(Ilow1:Ihigh1,Col.Lon), Cat(Ilow1:Ihigh1,Col.Lat));
        Dist2  = celestial.coo.sphere_dist_fast(Long,Lat, Cat(Ilow2:Ihigh2,Col.Lon), Cat(Ilow2:Ihigh2,Col.Lat));
        Ind    = [Ilow1-1+find(Dist1 <= Radius); Ilow2-1+find(Dist2 <= Radius)];
        
    else
        
        Inear = Util.find.mfind_bin(Cat(:,Col.Lon),[MinLong, MaxLong]);
        Ilow  = double(Inear(1));
        Ihigh = min(Ncat,double(Inear(2)) + 1);  % add 1 because of the way mfind_bin works
        Dist   = celestial.coo.sphere_dist_fast(Long,Lat, Cat(Ilow:Ihigh,Col.Lon), Cat(Ilow:Ihigh,Col.Lat));
        Ind    = Ilow-1+find(Dist <= Radius);
        
    end
end







