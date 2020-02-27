function [RA,Dec]=server_simbad(Name,OutUnits,CooType)
% Resolve an astronomical object name into coordinates using SIMBAD
% Package: VO.name
% Description: Resolve an astronomical object name into coordinates using
%              SIMBAD database.
% Input  : - String of object name (e.g., 'm81');
%          - Output units ['d','r','SH']. Default is 'd'.
%          - Coordinate type ['ICRS','FK5','FK4']. Default is 'ICRS'.
% Output : - J2000.0 RA [deg].
%          - J2000.0 Dec [deg].
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [RA,Dec]=VO.name.server_simbad('m31');
% Reliable: 2
%--------------------------------------------------------------------------

Def.OutUnits = 'd';
Def.CooType  = 'ICRS';
if (nargin==1)
    OutUnits = Def.OutUnits;
    CooType  = Def.CooType;
elseif (nargin==2)
    CooType  = Def.CooType;
elseif (nargin==3)
    % do nothing
else
    error('Illegal number of input arguments');
end

URL = sprintf('http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=%s&submit=SIMBAD+search',Name);

%Str = urlread(URL);
Str = webread(URL);

%Coo    = regexp(Str,'ICRS.*J2000.*TT*{1,15}(?<RA>\d\d\s\d\d\s\d\d\.\d\d\d)\s(?<Dec>[+-]\d\d\s\d\d\s\d\d\.\d\d).*<A HREF=','names');
%Coo    = regexp(Str,'ICRS.*J2000.*TT.{1,5}(?<RA>\d\d\s\d\d\s\d\d\.\d\d\d)\s(?<Dec>[+-]\d\d\s\d\d\s\d\d\.\d\d).*<A HREF=','names');
Coo    = regexp(Str,sprintf('%s.{1,300}TT.{1,5}(?<RA>\\d\\d\\s\\d\\d\\s\\d\\d\\.\\d{0,12})\\s(?<Dec>[+-]\\d\\d\\s\\d\\d\\s\\d\\d\\.\\d{0,11}).*<A HREF=',CooType),'names');
if (any(size(Coo) == 0))
    RA  = NaN;
    Dec = NaN;
else
    RA  = celestial.coo.convertdms(Coo.RA,'SHb',OutUnits);
    Dec = celestial.coo.convertdms(Coo.Dec,'SDb',OutUnits);
end
