function CatName=htmcat_names
% Get names of all HDF5/HTM catalogs in the /data/catsHTM/ directory.
% Package: VO.search
% Description: Get names of all HDF5/HTM catalogs in the /data/catsHTM/
%              directory.
% Input  : null
% Output : - A cell array of catalog names.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: CatName=VO.search.htmcat_names
% Reliable: 2
%--------------------------------------------------------------------------



P    = path;
RS   = regexp(P,':','split');
Icat = find(~Util.cell.isempty_cell(strfind(RS,'catsHTM')));
Ncat = numel(Icat);
CatName = cell(Ncat,1);
for I=1:1:Ncat
    Sp = regexp(RS{Icat(I)},filesep,'split');
    CatName{I} = Sp{end};
end