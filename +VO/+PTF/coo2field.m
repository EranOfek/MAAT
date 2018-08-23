function [Field,Offset]=coo2field(RA,Dec)
% Convert equatorial J2000 coordinates to PTF fields/CCDIDs
% Package: VO.PTF
% Description: Look for the PTF fields/CCDs that contains a given coordinate.
% Input  : - J2000.0 R.A. in radians or sexagesimal coordinates
%            (see convertdms.m for details).
%          - J2000.0 Dec. in radians or sexagesimal coordinates
%            (see convertdms.m for details).
% Output : - List of fields/CCDs that contains the coordinates:
%            [Field, CCD, RA, Dec, SDSS_coverage]
%            where RA and Dec are in degrees.
%          - Offset between coordinate and fields centers [RA, Dec]
%            in radians.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Jan 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Field,Offset]=VO.PTF.coo2field(1,1);
% Reliable: 2
%------------------------------------------------------------------------------

CatField = AstCat.CatField;
RAD = 180./pi;
Scale = VO.PTF.pixscale;
SizeX = 2048;
SizeY = 4096;

% Load the PTF Fields/CCD SDSS coverge into an AstCat object
Fields = AstCat.loadh2astcat('PTF_FieldCCD.hdf5');

RA  = celestial.coo.convertdms(RA,'gH','r');
Dec = celestial.coo.convertdms(Dec,'gD','R');

Col.RA  = colname2ind(Fields,'RA');
Col.Dec = colname2ind(Fields,'Dec');
%Res = search_cat(Fields.Cat(:,3:4),RA,Dec,'SearchRad',2./RAD,'CooType','sphere'); %,'ColX','RA','ColY','Dec');
Res = search_cat(Fields.(CatField)(:,[Col.RA, Col.Dec]),RA,Dec,'SearchRad',2./RAD,'CooType','sphere');

%D = sphere_dist(RA,Dec, Fields(:,Col.RA)./RAD, Fields(:,Col.Dec)./RAD);
% fields overlap candidates
%If = find(D<2./RAD);
If = Res.IndCat;
Nf = length(If);
[DelLon,DelLat]=celestial.coo.sphere_offset(RA.*ones(Nf,1),Dec.*ones(Nf,1),Fields.(CatField)(If,Col.RA),Fields.(CatField)(If,Col.Dec));

I = find(abs(DelLon)<0.5.*SizeX.*Scale./(RAD.*3600) & ...
         abs(DelLat)<0.5.*SizeY.*Scale./(RAD.*3600));


Field = Fields.(CatField)(If(I),:);
Offset = [DelLon(I), DelLat(I)];
