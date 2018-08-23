function [Vec,Res]=match_cats_pl(Cat,RefCat,varargin)
% Match two planer coordinates catalogs sorted by Y
% Package: VO.search
% Description: Given two planer coordinate catalogs. - for each entry
%              in the reference catalog (second input
%              argument), search for all nearby sources in the catalog
%              (first input).
% Input  : - A catalog sorted by Y.
%          - A reference catalog.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Radius' - Search radius. This is either a scalar or a vector
%                       which length is identical to that of the reference
%                       catalog (second input). If a vector than each
%                       source in the reference catalog may have a
%                       different search radius.
%                       Default is 2 (arcsec).
%            'CalcPA' - A logical flag indicating if to calculate and
%                       return also the position angle for each matched
%                       source.
%                       Default is false.
%            'ColRA'  - Index of logitude column in the first input
%                       catalog. Default is 1.
%            'ColDec' - Index of latitude column in the first input
%                       catalog. Default is 2.
%            'ColRefRA'-Index of logitude column in the second input
%                       reference catalog. Default is 1.
%            'ColRefDec'-Index of latitude column in the second input
%                       reference catalog. Default is 2.
%            'CheckIsSorted'- A logical flag indicating if to check that
%                       the first input catalog is sorted.
%                       Default is true.
% Output : - Structure that contains the following fields, each containing
%            a vector which length is identical to the size of the
%            reference catalog.
%            'Nfound' - Number of sources found in the catalog that within
%                       the search radius from the source in the reference
%                       catalog.
%            'MinDist'- Minimum distance (radians) of matched sources.
%                       NaN if not found.
%            'MinPA'  - P.A. of matched source with minimum distance.
%                       NaN if not found or 'CalcPA' is false.
%          - A structure array containing additional information on the
%            matched sources, an entry for each matched refernce source,
%            with the following fields:
%            'IndRef' - Index of source in reference catalog.
%            'IndCat' - List of indices in the catalog that are matched to
%                       the 'IndRef' source in the reference catalog.
%            'Dist'   - Vecor of angular distances (radians) for each one
%                       of the sources indicated in 'IndCat'.
%            'PA'     - Vector of position angles (radians).
%            'Num'    - Number of sources within search radius.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=sortrows(rand(1000,2),2); Ref=rand(1000,2);
%          Cat = Cat.*5000; Ref = Ref.*5000;
%          Cat(1:50,:) = Ref(1:50,:);
%          Cat = sortrows(Cat,2);
%          [Vec,Res]=VO.search.match_cats_pl(Cat,Ref)
% Reliable: 
%--------------------------------------------------------------------------


DefV.Radius               = 2;
DefV.CalcPA               = false;
DefV.ColRA                = 1;
DefV.ColDec               = 2;
DefV.ColRefRA             = 1;
DefV.ColRefDec            = 2;
DefV.CheckIsSorted        = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (InPar.CheckIsSorted)
    if (~issorted(Cat(:,InPar.ColDec)))
        error('Input catalog must be sorted by Y axis');
    end
end


% number of sources
Ncat = size(Cat,1);
Nref = size(RefCat,1);

Radius = InPar.Radius(:).'.*ones(1,Nref);

% for each source in the reference
Res = Util.struct.struct_def({'IndRef','IndCat','Dist','PA','Num'},0,1);
    
% binary search on sorted catalog
% search simultanously all ref sources
Iupp = Util.find.mfind_bin(Cat(:,InPar.ColDec),RefCat(:,InPar.ColRefDec).'+Radius);
Ilow = Util.find.mfind_bin(Cat(:,InPar.ColDec),RefCat(:,InPar.ColRefDec).'-Radius);


Ilow = max(1,double(Ilow));
Iupp = min(Ncat,double(Iupp + 1));  % add 1 because of the way mfind_bin works

Ncand = Iupp-Ilow;
Ic    = find(Ncand>=1);
Nc    = numel(Ic);

Vec.Nfound  = zeros(Nref,1);
Vec.MinDist = nan(Nref,1);
Vec.MinPA   = nan(Nref,1);
K = 0;

for Icr=1:1:Nc
    Iref = Ic(Icr);  % single
    Icat = (Ilow(Iref):Iupp(Iref)); % multiple

    if (InPar.CalcPA)
        [Dist,PA]  = Util.Geom.plane_dist(Cat(Icat,InPar.ColRA),...
                                             Cat(Icat,InPar.ColDec),...
                                             RefCat(Iref,InPar.ColRefRA),...
                                             RefCat(Iref,InPar.ColRefDec));
    else
        
        [Dist]  = Util.Geom.plane_dist(Cat(Icat,InPar.ColRA),...
                                             Cat(Icat,InPar.ColDec),...
                                             RefCat(Iref,InPar.ColRefRA),...
                                             RefCat(Iref,InPar.ColRefDec));
       
    end
                                         
    IndRelative = find(Dist <= Radius(Iref));              
    IndCat      = Ilow(Icr)-1+IndRelative;

    Vec.Nfound(Iref) = numel(IndCat);
    if (Vec.Nfound(Iref)>0)
        
        [Vec.MinDist(Iref),MinInd] = min(Dist);
        if (InPar.CalcPA)
            Vec.MinPA(Iref)            = PA(MinInd);
        end
        
        % set the additonal output info
        if (nargout>1)
            K = K + 1;
            Res(K).IndCat = IndCat;
            Res(K).IndRef = Iref;
            Res(K).Num    = numel(IndCat);
            Res(K).Dist   = Dist(IndRelative);
            
            if (InPar.CalcPA)
                [~,Res(K).PA]  = Util.Geom.plane_dist(Cat(Iref,InPar.ColRA),...
                                                 Cat(Iref,InPar.ColDec),...
                                                 RefCat(IndCat,InPar.ColRefRA),...
                                                 RefCat(IndCat,InPar.ColRefDec));
            end
        end
    end
        
end

    
     