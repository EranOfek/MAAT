function Link=navigator_link(RA,Dec,Size)
% Given J2000 equatorial coordinates get link to PS1 navigator image.
% Package: telescope.PS1
% Description: Get link to PS1 navigator image
% Input  : - J2000.0 R.A. [rad, [H M S], sexagesimal string]
%          - J2000.0 Dec. [rad, [Sign D M S], sexagesimal string]
%          - Cutout size [arcsec]. Default is 240.
% Output : - Cell array of links.
%     By : Eran O. Ofek                    Dec 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Link=telescope.PS1.navigator_link(pi,0.4)
% Reliable: 2

if (nargin<3)
    Size = 240;
end

RA  = celestial.coo.convertdms(RA,'gH','SH');
Dec = celestial.coo.convertdms(Dec,'gD','SD');

N = size(RA,1);
Link = cell(N,1);
for I=1:1:N
    Link{I} = sprintf('http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=%s%s&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size=%d&output_size=0&verbose=0&autoscale=99.500000&catlist=',...
                   RA(I,:),Dec(I,:),Size);
end

