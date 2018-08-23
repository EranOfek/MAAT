function Link=ned_link(RA,Dec,Radius)
% Get link to NED search by coordinates
% Package: VO.NED
% Description: Get link to NED search by coordinates
% Input  : - RA [rad] or sexagesimal.
%          - Dec [rad] or sexagesimal
%          - Radius [arcmin]. Default is 3.
% Output : - Cell array of NED links.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

if (nargin<3)
    Radius = 3; % [arcmin]
end

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

RAv  = celestial.coo.convertdms(RA,'gH','SH');
Decv = celestial.coo.convertdms(Dec,'gD','SD');
Ncoo = size(RAv,1);

Link = cell(Ncoo,1);
for Icoo=1:1:Ncoo
    Link{Icoo} = sprintf('https://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon=%s&lat=%s&radius=%f&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&search_type=Near+Position+Search&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES',...
                    RAv(Icoo,:),...
                    Decv(Icoo,:),...
                    Radius);
end

