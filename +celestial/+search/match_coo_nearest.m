function [Ind,Nmatch,Dist,PA,FlagMatched2]=match_coo_nearest(Data1,Data2,Radius,varargin)
% Match two lists by spherical coordinates for nearest source only.
% Package: celestial.search
% Description: Given two lists with spherical coordinates, search for
%              nearest object in the second list associated with each
%              object in the first list.
% Input  : - First list. By default [Lon, Lat]. Coordinates in radians.
%          - Second list sorted by lat. By default [Lon, Lat].
%            Coordinates in radians.
%          - Match radius [radians]. Scalar or vector.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MinDist'     - Minimum distance [radians].
%                            Objects with distance smaller than this
%                            distance will be excluded.
%                            Default is 0.
%            'CheckSorted' - Verify that second list is sorted by Lat.
%                            Default is false.
%            'ColRA1'      - Column of longitude in 1st list.
%            'ColDec1'     - Column of latitude in 1st list.
%            'ColRA2'      - Column of longitude in 2nd list.
%            'ColDec2'     - Column of latitude in 2nd list.
% Output : - A vecor of indices of the nearest source in the 2nd list found
%            within the angular distance from a source in the 1st list.
%          - Vector of the number of matched per source.
%          - Like 1st output, but for angular distance [radians].
%          - Like 1st output, but for position angle [radians].
%          - Flag indicating if an entry in the 2nd list was matched.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=celestial.search.match_coo_nearest(Data1,Data2,0.001);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.MinDist              = 0;
DefV.CheckSorted          = false;
DefV.ColRA1               = 1;
DefV.ColDec1              = 2;
DefV.ColRA2               = 1;
DefV.ColDec2              = 2;
if (isempty(varargin))
    InPar = DefV;
else
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
end

if (InPar.CheckSorted)
    if (~issorted(Data2(:,InPar.ColDec2)))
        error('Data should be sorted by Dec');
    end
end

Im1 = Util.find.mfind_bin(Data2(:,InPar.ColDec2),Data1(:,InPar.ColDec1).'-Radius.');
Im2 = Util.find.mfind_bin(Data2(:,InPar.ColDec2),Data1(:,InPar.ColDec1).'+Radius.');
[N2,Ncol2]  = size(Data2);
I1  = max(Im1-1,1);
I2  = min(Im2+1,N2);

[N1,Ncol1]  = size(Data1);

if (numel(Radius)==1)
    Radius = Radius.*ones(N1,1);
end

Ind    = nan(N1,1);
Nmatch = zeros(N1,1);
Dist   = nan(N1,1);
PA     = nan(N1,1);


for Id1=1:1:N1
    
    D    = celestial.coo.sphere_dist_fast(Data2(I1(Id1):I2(Id1),InPar.ColRA2),Data2(I1(Id1):I2(Id1),InPar.ColDec2),...
                                          Data1(Id1,InPar.ColRA1),Data1(Id1,InPar.ColDec1));

    II  = uint32(find(D<Radius(Id1) & D>InPar.MinDist));
    if (~isempty(II))
        [MinD,MinI] = min(D(II));
        Ind(Id1)    = I1(Id1) + II(MinI) - 1;
        Nmatch(Id1) = numel(II);
        
        if (nargout>2)
            Dist(Id1) = MinD;
            if (nargout>3)
                [~,PA(Id1)] = celestial.coo.sphere_dist_fast(Data2(Ind(Id1),InPar.ColRA2),Data2(Ind(Id1),InPar.ColDec2),...
                                              Data1(Id1,InPar.ColRA1),Data1(Id1,InPar.ColDec1));
            end
        end
    end
end

if (nargout>3)
    Matched2 = Ind;
    Matched2 = Matched2(Matched2~=0);
    Matched2 = unique(Matched2);
    FlagMatched2 = ismember((1:1:N2),Matched2);
end


