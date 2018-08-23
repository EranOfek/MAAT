function [LinkCutUnm,LinkFull,Link]=get_stack(RA,Dec,Size)
% Get links to PS1 corrected stack FITS images.
% Package: telescope.PS1
% Description: Get link to PS1 FITS images.
% Input  : - J2000.0 R.A. [rad, [H M S], sexagesimal string]
%          - J2000.0 Dec. [rad, [Sign D M S], sexagesimal string]
%          - Cutout size [arcsec]. Default is 240.
% Output : - A 5 column cell array of links to the cutout images.
%            Each column for specific filter (g,r,i,z,y) and each line
%            for specific requested coordinate.
%          - A 5 column cell array of links to the full images.
%            Each column for specific filter (g,r,i,z,y) and each line
%            for specific requested coordinate.
%          - Cell array of links to the navigator pages.
%     By : Eran O. Ofek                    Dec 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [LinkCut,LinkFull,Link]=VO.PS1.get_stack(pi,0.4)
% Reliable: 2

BaseURL = 'http://ps1images.stsci.edu/rings';

if (nargin<3)
    Size = 240;
end

Ncolor = 5;

Link=VO.PS1.navigator_link(RA,Dec,Size);
Nlink = numel(Link);

LinkCut  = cell(Nlink,Ncolor);
LinkFull = cell(Nlink,Ncolor);


for Ilink=1:1:Nlink
    List = www.find_urls(Link{Ilink},'match','.*?\.fits$');
    
    % full images
    IndFull = ~Util.cell.isempty_cell(strfind(List,BaseURL));
    [LinkFull{Ilink,:}] = deal(List{IndFull});
    
    % cut images
    %IndCut = ~Util.cell.isempty_cell(strfind(List,'http://ps1images.stsci.edu/cgi-bin/fitscut'));
    IndCut = ~Util.cell.isempty_cell(strfind(List,'cutout_rings'));
    [LinkCut{Ilink,:}] = deal(List{IndCut});
    
end
LinkCut = regexprep(LinkCut,'//ps1images.stsci.edu/cgi-bin/','/cgi-bin/');

LinkCutUnm = LinkCut;



Tmp = regexp(LinkCut,'&amp;','split');
N = numel(Tmp);
for I=1:1:N
    LinkCut{I} = Tmp{I}{1};
end
    
LinkCut = regexprep(LinkCut,'http://ps1images.stsci.edu//','http://');

