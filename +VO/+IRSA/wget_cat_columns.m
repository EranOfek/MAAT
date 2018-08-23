function [Cols]=wget_cat_columns(CatName)
% Retrieve a list of all columns for an IPAC/IRSA public catalog
% Package: VO.IRSA
% Description: Given an IPAC/IRSA catalog name
%             (see Catalog.IRSA.wget_all_catnames) retrieve all columns
%             in the catalog.
% Input  : - Catalog name (e.g., 'wise_allwise_p3as_psd').
%            Use Catalog.IRSA.wget_all_catname to retrieve a list of
%            catalogs.
% Output : - A structure array containing the column names and information.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: http://irsa.ipac.caltech.edu/applications/Gator/GatorAid/irsa/catsearch.html
% Example: [Cols]=VO.IRSA.wget_cat_columns('wise_allwise_p3as_psd');
% Reliable: 2
%--------------------------------------------------------------------------

Url = VO.IRSA.irsa_db_url;

Tbl = urlread(sprintf('%snph-dd?mode=ascii&catalog=%s',Url,CatName));
TblCat = VO.IRSA.read_ipac_table(Tbl,'str');  

Cols.Name   = Util.string.spacedel(TblCat.Cat(:,TblCat.Col.name));
Cols.Desc   = strtrim(TblCat.Cat(:,TblCat.Col.description));
Cols.Units  = strtrim(TblCat.Cat(:,TblCat.Col.units));
Cols.Format = strtrim(TblCat.Cat(:,TblCat.Col.format));
Cols.Nulls  = strtrim(TblCat.Cat(:,TblCat.Col.nulls));
