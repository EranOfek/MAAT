function Link=skymapper_cutout_link(RA,Dec,varargin)
% Generate a URL link to SkyMapper image cutouts.
% Package: VO.SkyMapper
% Description: Generate a URL link to SkyMapper image cutouts.
% Input  : - J2000.0 R.A. [rad] or see celestial.coo.convertdms.
%          - J2000.0 Dec. [rad] or see celestial.coo.convertdms.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Size' - Cutout size [deg deg]. Default is [0.1 0.1].
%            'Format' - 'fits' | 'png'. Default is 'fits'.
%            'Bands'  - Defaukt is 'y,v,g,r,i,z'.
%            'Intersect - 'covers' | 'center' | 'overlaps'
%                       Default s 'covers'.
%            'StartJD' - Default is [1 1 2000].
%            'EndJD'   - Default is [1 1 2100].
%            'BaseURL' - Default is
%                   'http://skymappersiap.asvo.nci.org.au/dr1_cutout/query?'.
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Link=VO.SkyMapper.skymapper_cutout_link(1,-1)
% Reliable: 2
%--------------------------------------------------------------------------

DefV.Size                 = [0.1 0.1];  % [deg]
DefV.Format               = 'fits';     % fits | png
DefV.Bands                = 'y,v,g,r,i,z';
DefV.Intersect            = 'covers';   % covers | center | overlaps
DefV.StartJD              = [1 1 2000];
DefV.EndJD                = [1 1 2100];
DefV.BaseURL              = 'http://skymappersiap.asvo.nci.org.au/dr1_cutout/query?';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (numel(InPar.StartJD)==1)
    StartMJD = celestial.time.jd2mjd(InPar.StartJD);
else
    StartMJD = celestial.time.jd2mjd(celestial.time.julday(InPar.StartJD));
end
if (numel(InPar.EndJD)==1)
    EndMJD = celestial.time.jd2mjd(InPar.EndJD);
else
    EndMJD = celestial.time.jd2mjd(celestial.time.julday(InPar.EndJD));
end

RA  = celestial.coo.convertdms(RA,'gH','d');
Dec = celestial.coo.convertdms(Dec,'gD','d');


SizeStr   = sprintf('SIZE=%f,%f',InPar.Size.*ones(1,2));
FormatStr = sprintf('FORMAT=image/%s',InPar.Format);
BandStr   = sprintf('BAND=%s',InPar.Bands);
InterStr  = sprintf('INTERSECT=%s',InPar.Intersect);
StartStr  = sprintf('MJD_START=%f',StartMJD);
EndStr    = sprintf('MJD_END=%f',EndMJD);

Ncoo = numel(RA);
Link = cell(Ncoo,1);
for Icoo=1:1:Ncoo
    PosStr     = sprintf('POS=%f,%f',RA(Icoo),Dec(Icoo));
    Link{Icoo} = sprintf('%s%s&%s&%s&%s&%s&%s&%s',InPar.BaseURL,...
                    PosStr,SizeStr,...
                    FormatStr, BandStr, InterStr,...
                    StartStr, EndStr);
end
