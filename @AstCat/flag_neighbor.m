function FlagNeigh=flag_neighbor(AstC,varargin)
% Flag sources with neighbors above some flux ratio threshold
% Package: @AstCat
% Package Category: search
% Description: Given an AstCat object with a catalog containing X and Y
%              coordinates and flux ratio flag sources which have neighbors
%              with brighter than some flux ratio relative to the source.
%              The user can define a list of distance and flux ratio
%              thresholds.
% Input  : - An AstCat object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'NeighRadius' -  Vector of radii to search. Each radius is
%                             associated with a flux ratio threshold.
%                             Default is [10 30].
%            'NeighMaxFluxRatio' - Vector of flux ratio threshold
%                             corresponding to each radius.
%                             Default is [0.1 4].
%            'FluxColName' - Flux column name in the AstCat object. If
%                            empty then will use only radius search.
%                            Default is 'FLUX_PSF'.
%            'XColName'    - X coordinate column name.
%                            Default is 'XWIN_IMAGE.
%            'YColName'    - Y coordinate column name.
%                            Default is 'YWIN_IMAGE.
% Output : - A vector of logicals indicating if each source have a
%            neighbor.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FlagNeigh=flag_neighbor(Sim)
% Reliable: 2

DefV.NeighRadius         = [5 20]; %[10 30];
DefV.NeighMaxFluxRatio   = [0.3 2];
DefV.FluxColName         = 'SN';
DefV.XColName            = 'XWIN_IMAGE';
DefV.YColName            = 'YWIN_IMAGE';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Ncat = numel(AstC);
if (Ncat>1)
    error('Works on a single element AstCat object');
end

Nradius   = numel(InPar.NeighRadius);
MaxRadius = max(InPar.NeighRadius);

IsSorted = issorted(AstC,InPar.YColName);
if (~all(IsSorted))
    error('AstCat ocatalogs must be sorted by Y');
end

for Icat=1:1:Ncat
    % for each catalog
    
    X    = col_get(AstC(Icat),InPar.XColName);
    Y    = col_get(AstC(Icat),InPar.YColName);
    if (~isempty(InPar.FluxColName))
        Flux = col_get(AstC(Icat),InPar.FluxColName);
    end
    Nsrc = numel(X);
    IndLow  = Util.find.mfind_bin(Y,Y'-MaxRadius);
    IndHigh = Util.find.mfind_bin(Y,Y'+MaxRadius);
    
    FlagNeigh = false(Nsrc,1);
    for Isrc=1:1:Nsrc
        IndRange  = setdiff((IndLow(Isrc):IndHigh(Isrc)),Isrc);
        Dist      = sqrt((X(Isrc)-X(IndRange)).^2 + (Y(Isrc)-Y(IndRange)).^2);
        if (isempty(InPar.FluxColName))
            for Iradius=1:1:Nradius
                FlagNeigh(Isrc) = FlagNeigh(Isrc) || any(Dist<InPar.NeighRadius(Iradius));
            end
        else
            FluxRatio = Flux(Isrc)./Flux(IndRange);

            for Iradius=1:1:Nradius
                FlagNeigh(Isrc) = FlagNeigh(Isrc) || any(Dist<InPar.NeighRadius(Iradius) & FluxRatio<InPar.NeighMaxFluxRatio(Iradius));
            end
        end
    end
end
            
            
    