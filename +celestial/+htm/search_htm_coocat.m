function HTMind=search_htm_coocat(RA,Dec,SearchRad,HTMcenters,Use_search_cat)
%--------------------------------------------------------------------------
% search_htm_coocat function                                           htm
% Description: 
% Input  : - RA (J2000.0) (radians, sexahesimal string or [H M S]).
%          - Dec (J2000.0) (radians, sexahesimal string or [Sign D M S]).
%          - Search radius [radians].
%          - Either a file name containing the catalog of HTM centers or
%            the HTM centers.
%            The HTM centers is a structure array with a 'Cat' field
%            containing a 3 columns matrix of [RA, Dec, Ptr], sorted by
%            the declination
%          - Use search_cat.m (true) or use simple sphere_dist (false.
%            Use search_cat.m for large catalogs. Default is true.
% Output : - Indices of the low level HTM in which the search coordinates
%            may reside in, and there are more than 0 sources.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: HTMind=celestial.htm.search_htm_coocat(1,1,1./RAD,'FIRST_htm.mat')
% Reliable: 2
%--------------------------------------------------------------------------
import Util.IO.*

InvRAD    = pi./180;

CatField  = 'Cat';
ColField  = 'Col';

if (nargin<5)
    Use_search_cat = true;
end


if (ischar(HTMcenters))
    HTMcenters = load2(HTMcenters);
end

ColRA     = HTMcenters.(ColField).RA;
ColDec    = HTMcenters.(ColField).Dec;
ColPtr    = HTMcenters.(ColField).Ptr;
if (isfield(HTMcenters.(ColField),'Nsrc'))
    ColNsrc   = HTMcenters.(ColField).Nsrc;
else
    ColNsrc   = [];
end
HTMcenters.(CatField) = double(HTMcenters.(CatField));

Ncen = size(HTMcenters.(CatField),1);  % number of HTM in lower level
Nlev = round(log(Ncen.*0.5)./log(4));  % number of levels
TriSide = 180./(2.^Nlev);              % estimated size of the triangle side
TriCM   = TriSide.*(sqrt(2)./2).*1.01;         % estimated size of distance from tri edge to center of mass
TriCM   = TriCM.*InvRAD;               % [radians]

if (Use_search_cat)
    Res     = search_cat(HTMcenters.(CatField)(:,[ColRA,ColDec]),RA,Dec,...
                         'SearchRad',SearchRad+TriCM,'SearchMethod','binms','IsRad',true);
    HTMind  = HTMcenters.(CatField)(Res.IndCat,ColPtr);
else
    D = celestial.coo.sphere_dist_fast(HTMcenters.(CatField)(:,ColRA),HTMcenters.(CatField)(:,ColDec),...
                         RA,Dec);
    if (isempty(ColNsrc))
        % do not check if Nsrc>0
        HTMind  = HTMcenters.(CatField)(D<(SearchRad+TriCM),ColPtr);
    else
        % check if Nsrc>0
        HTMind  = HTMcenters.(CatField)(D<(SearchRad+TriCM) & HTMcenters.(CatField)(:,ColNsrc)>0 ,ColPtr);
    end
end

