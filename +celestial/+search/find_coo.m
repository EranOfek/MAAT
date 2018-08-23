function [Ind,DistI,PAI]=find_coo(Data,RA,Dec,Radius,varargin)
% Cone search in a table with spherical coordinates.
% Package: celestial.search
% Description: Search for a coordinate within a radius in a table of
%              spherical coordinates.
% Input  : - Table sorted by latitude. By default [Long, Lat] in radians.
%          - Longitude to search [radians] or sexagesimal string.
%          - Latitude to search [radians] or sexagesimal string.
%          - Search radius [radians].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'CheckSorted' - Check if table is sorted. Default is false.
%            'ColRA'       - Column of longitude. Default is 1.
%            'ColDec'      - Column of latitude. Default is 1.
% Output : - Indices of entries found in search radius.
%          - Distance [radians] of found entries from search point.
%          - Position angle [radians] of found entries from search point.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Ind,DistI,PAI]=celestial.search.find_coo(rand(10,2),0.4,0.4,0.1)
% Reliable: 2
%--------------------------------------------------------------------------

DefV.CheckSorted          = false;
DefV.ColRA                = 1;
DefV.ColDec               = 2;
if (isempty(varargin))
    InPar = DefV;
else
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
end

if (InPar.CheckSorted)
    if (~issorted(Data(:,InPar.ColDec)))
        error('Data should be sorted by Dec');
    end
end

if (ischar(RA))
    RA = celestial.coo.convertdms(RA,'SH','r');
end

if (ischar(Dec))
    Dec = celestial.coo.convertdms(Dec,'SD','R');
end

I1 = Util.find.bin_sear(Data(:,InPar.ColDec),Dec-Radius);
I2 = Util.find.bin_sear(Data(:,InPar.ColDec),Dec+Radius);
N  = size(Data,1);
I1 = max(I1-1,1);
I2 = min(I2+1,N);

D   = celestial.coo.sphere_dist_fast(Data(I1:I2,InPar.ColRA),Data(I1:I2,InPar.ColDec),RA,Dec);
II  = find(D<Radius);
Ind = I1 + II - 1;

if (nargout>1)
    DistI = D(II);
    
    if (nargout>2)
        [~,PAI] = celestial.coo.sphere_dist_fast(Data(Ind,InPar.ColRA),Data(Ind,InPar.ColDec),RA,Dec);
    end
end



