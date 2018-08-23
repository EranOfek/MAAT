function [Data,ColCell,C,Flag]=gaia_dr1_read_file(FileName)
% Read GAIA-DR1 file for reformatting purposses
% Description: Given a GAIA-DR1 CSV file name, read the file into a cell
%              array and also return a matrix of selected columns of
%              good entries.
% Input  : - File name.
% Output : - Matrix of selected columns of all entries
%          - Cell array of column names:
%            [RA, ErrRA, Dec, ErrDec, MagG, ErrMagG].
%          - Cell array of all columns.
%          - Flag of good entries
% Example: cd /raid/eran/Catalogue/GAIA-DR1/Orig
%          Data=VO.prep..GAIA.gaia_dr1_read_file('GaiaSource_000-012-071.csv');

[ColCell,Col]=VO.prep.GAIA.gaia_dr1_cat_columns;
Ncol   = numel(ColCell);
Format1 = str_duplicate('%f ',33,'');
Format2 = str_duplicate('%f ',5,'');
Format3 = str_duplicate('%f ',12,'');
Format4 = str_duplicate('%f ',4,'');


Format = sprintf('%s %%s %s %%s %s %%s %s\n',Format1,Format2,Format3,Format4);

FID = fopen(FileName,'r');
C = textscan(FID,Format,'Delimiter',',','HeaderLines',1);
fclose(FID);


% selection criteria:
% astrometric_excess_noise (32) < 5 mas
% ra_err (6) < 5 mas
% dec_err (8) < 5 mas
% C{51}./C{50} < 0.02  % rel phot err
% ra_dec_corr (15) -0.95..0.95

ColCell = {'RA','ErrRA','Dec','ErrDec','MagG','ErrMagG','ExcessNoise','ExcessNoiseSig'};
Flag = C{6}<5 & C{8}<5 & C{32}<5 & C{51}./C{50}<0.02 & C{15}>-0.95 & C{15}<0.95;
Data = double([C{5}, C{6}, C{7}, C{8}, C{52}, 1.086.*C{51}./C{50}, C{32}, C{33}]);


