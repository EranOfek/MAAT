function [Cat,ColCell,Col]=wget_sdss(RA,Dec,SearchRad,Shape,ColCell,CooOut)
% Query SDSS PhotoPrimary table around specific coordinate.
% Package: VO.SDSS
% Description: Query SDSS PhotoPrimary table around specific coordinate.
%              See VO.SDSS.run_sdss_sql.m for a more general queries.
% Input  : - J2000.0 RA (rad, sexagesimal string or [H M S]).
%          - J2000.0 Dec (rad, sexagesimal string or [Sign D M S]).
%          - Search radius (radians). Default is 1 deg.
%          - Search shape {'circ'|'box'}. Default is 'circ';
%          - Cell array of PhotoPrimary columns to retrieve.
%            Default is:
%            {'RA','Dec','type','modelMag_u','modelMagErr_u','modelMag_g',
%             'modelMagErr_g','modelMag_r','modelMagErr_r','modelMag_i',
%             'modelMagErr_i','modelMag_z','modelMagErr_z'};
%            If empty then use default.
%          - Units of output ra and dec {'rad','deg'}. Default is 'rad'.
% Output : - Matrix of SDSS catalog.
%          - Cell array of column names.
%          - Structure array of column indices.
% See also: run_sdss_sql.m
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cat,ColCell,Col]=VO.SDSS.wget_sdss('12:00:00','+30:00:00',0.1./RAD);
% Reliable: 2
%--------------------------------------------------------------------------
RAD    = 180./pi;
InvRAD = pi./180;


Def.SearchRad = InvRAD;
Def.Shape     = 'circ';
Def.ColCell   = {'ra','dec','type','flags','modelMag_u','modelMagErr_u','modelMag_g','modelMagErr_g','modelMag_r','modelMagErr_r','modelMag_i','modelMagErr_i','modelMag_z','modelMagErr_z'};
Def.CooOut    = 'rad';

if (nargin==2)
    SearchRad   = Def.SearchRad;
    Shape       = Def.Shape;
    ColCell     = Def.ColCell;
    CooOut      = Def.CooOut;
elseif (nargin==3)
    Shape       = Def.Shape;
    ColCell     = Def.ColCell;
    CooOut      = Def.CooOut;
elseif (nargin==4)
    ColCell     = Def.ColCell;
    CooOut      = Def.CooOut;
elseif (nargin==5)
    CooOut      = Def.CooOut;
elseif (nargin==6)
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isempty(ColCell))
    ColCell = Def.ColCell;
end


RA  = celestial.coo.convertdms(RA,'gH','r');
Dec = celestial.coo.convertdms(Dec,'gD','R');


if (numel(SearchRad)==1)
    SearchRad = [SearchRad, SearchRad];
end

Q{1} = ColCell;
Q{2} = {'PhotoPrimary'};

MinDec = (Dec - SearchRad(2)).*RAD;
MaxDec = (Dec + SearchRad(2)).*RAD;

MinRA  = (RA  - SearchRad(1)./cosd(MaxDec)).*RAD;
MaxRA  = (RA  + SearchRad(1)./cosd(MaxDec)).*RAD;

if (MinRA<0)
    Q{3} = sprintf('((ra between %10.6f and %10.6f) or (ra between %10.6f and %10.6f)) and (dec between %10.6f and %10.6f)',...
                   360+MinRA,360,0,MaxRA,MinDec,MaxDec);
else
    if (MaxRA>360)
        Q{3} = sprintf('((ra between %10.6f and %10.6f) or (ra between %10.6f and %10.6f)) and (dec between %10.6f and %10.6f)',...
                       MinRA,360,0,MaxRA-360,MinDec,MaxDec);
    else
        Q{3} = sprintf('(ra between %10.6f and %10.6f) and (dec between %10.6f and %10.6f)',...
                       MinRA,MaxRA,MinDec,MaxDec);
    end
end


[Cat,Msg] = VO.SDSS.run_sdss_sql(Q);
if (~isempty(Msg))
    Msg
    error('run_sdss_sql.m may have failed');
end
Col = cell2struct(num2cell(1:1:length(ColCell)),ColCell,2);
Col.RA  = Col.ra;    % synonyms
Col.Dec = Col.dec;

switch lower(Shape)
    case {'circle','circ'}
        if (~isempty(Cat))
            D = celestial.coo.sphere_dist_fast(Cat(:,Col.ra).*InvRAD,Cat(:,Col.dec).*InvRAD,RA,Dec);
            Cat = Cat(D<=SearchRad(1),:);
        end
    otherwise
        % assume box - do nothing
end

switch lower(CooOut)
    case 'rad'
        Cat(:,[Col.ra, Col.dec]) = Cat(:,[Col.ra, Col.dec]).*InvRAD;
    otherwise
        % assume requested output is 'deg'
        % do nothing
end

% sort by declination
Cat = sortrows(Cat,Col.dec);
        