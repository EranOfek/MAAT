function [AllCat,Desc]=wget_all_catnames
% Retrieve a list of all IPAC/IRSA public catalogs
% Package: VO.IRSA
% Description: Retrieve a list of all IPAC/IRSA public catalogs
% Input  : null
% Output : - A cell array of all public catalog names.
%          - A cell array of catalog descriptions.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: http://irsa.ipac.caltech.edu/applications/Gator/GatorAid/irsa/catsearch.html
% Example: [AllCat,Desc]=VO.IRSA.wget_all_catnames; [AllCat,Desc]
% Reliable: 2
%--------------------------------------------------------------------------

Url = VO.IRSA.irsa_db_url;
Tbl = urlread(sprintf('%snph-scan?mode=ascii',Url));

TblCat  = VO.IRSA.read_ipac_table(Tbl,'str');
CatName = TblCat.Cat(:,TblCat.Col.catname);
AllCat  = Util.string.spacedel(CatName);

Desc    = strtrim(TblCat.Cat(:,TblCat.Col.description));

