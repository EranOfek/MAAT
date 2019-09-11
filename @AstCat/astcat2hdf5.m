function astcat2hdf5(AC,FileName,varargin)
% Save AstCat object in HDF5 file
% Package: @AstCat
% Description: Save AstCat object in HDF5 file.
%              If the AstCat object has multiple elements, then each table
%              is saved as a different dataset.
% Input  : - An AstCat object.
%          - HDF5 file name to create.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'DataName' - A cell array of data names.
%                         If empty, then use default (i.e., V1, V2,...).
%                         Default is empty.
%            'Convert2Real' - Default is true.
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

CatField = AstCat.CatField;



DefV.DataName             = []; 
DefV.Convert2Real         = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

N = numel(AC);

for I=1:1:N
    if isempty(InPar.DataName)
        DN = sprintf('/V%d',I);
    else
        DN = InPar.DataName{I};
    end
    ColCell = [AC(I).ColCell(:),num2cell(1:1:numel(AC(I).ColCell))'];
    if (InPar.Convert2Real)
        AC(I).(CatField) = real(AC(I).(CatField));
    end
    HDF5.save(AC(I).(CatField),FileName,DN,ColCell);
end