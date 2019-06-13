function [RA,Dec]=server_ned(Name,OutUnits)
% Resolve an astronomical object name into coordinates using NED.
% Package: VO.name
% Description: Resolve an astronomical object name into coordinates using
%              NASA Extragalactic Database (NED).
% Input  : - String of object name (e.g., 'm81');
%          - Output units {'d','r','SH'}. Default is 'd'.
% Output : - J2000.0 RA [deg].
%          - J2000.0 Dec [deg].
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [RA,Dec]=VO.name.server_ned('m31');
% Reliable: 2
%--------------------------------------------------------------------------

Def.OutUnits = 'd';
if (nargin==1)
    OutUnits = Def.OutUnits;
end

URL = sprintf('http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES',Name);

%Str = urlread(URL);
Str = webread(URL);

%SubStr = regexp(Str,'<A HREF="#ObjNo1">1.*\d\dh\d\dm\d\d\.\ds\s[+-]\d\dd\d\dm\d\ds.{1,100}(<A HREF)','match');
Coo    = regexp(Str,'<A HREF="#ObjNo1">1.*(?<RA>\d\dh\d\dm\d\d\.\ds)\s(?<Dec>[+-]\d\dd\d\dm\d\ds).{1,100}(<A HREF)','names');
if (~isempty(Coo))
    RA  = celestial.coo.convertdms(Coo.RA,'SHh',OutUnits);
    Dec = celestial.coo.convertdms(Coo.Dec,'SDh',OutUnits);
else
    RA  = NaN;
    Dec = NaN;
end
