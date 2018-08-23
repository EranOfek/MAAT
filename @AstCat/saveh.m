function saveh(AstC,FileName)
%--------------------------------------------------------------------------
% saveh function                                             class/@AstCat
% Description: Save an AstCat object in an HDF5 file.
%              The Column names and indices will be written as attributes.
% Input  : - An AstCat object.
%          - HDF5 file name in which to write the catalogs.
% Output : null.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Aug 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: saveh(AstC,'file.hdf5')
% Reliable: 2
%--------------------------------------------------------------------------

CatField = AstCat.CatField;
ColField = AstCat.ColField;

Ncat = numel(AstC);
for Icat=1:1:Ncat
    % for each catalog
    %if (istable(AstC(Icat).(CatField))),
    %    error('No support to table format yet');
    %end
    VarName = sprintf('Cat%d',Icat);
    
    AttribCell = [fieldnames(AstC(Icat).Col), struct2cell(AstC(Icat).(ColField))].';
    Util.IO.saveh(FileName,AstC(Icat).(CatField),VarName,AttribCell); % Na'ama, 20180807
end