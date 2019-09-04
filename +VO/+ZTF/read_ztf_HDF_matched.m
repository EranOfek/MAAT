function [Cat,ColCell,CatProp]=read_ztf_HDF_matched(FieldID,Lines,varargin)
% Read ZTF matched light curves from local HDF5 light curve files.
% Package: VO.ZTF
% Description: Read ZTF matched light curves from local HDF5 light curve
%              files. The HDF5 files are distributed as part of the catsHTM
%              catalogs.
% Input  : - ZTF field number.
%          - [start end] lines to read. The lines for a given source are
%            available in I1 and I2 in the 'ztfSrcLCDR1' catsHTM catalog.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FileTemp' - File template name. Default is
%                       'ztfLCDR1_%06d.hdf5'.
%            'ColCell'  - Column names for catalog.
%                       Default is
%                       {'HMJD','Mag','MagErr','ColorCoef','Flags'}.
% Output : - Catalog
%          - ColCell
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=VO.ZTF.read_ztf_HDF_matched(686,[1 100])
%          Cat=VO.ZTF.read_ztf_HDF_matched(703,[38104798    38104901])
% Reliable: 2
%--------------------------------------------------------------------------


DefV.FileTemp           = 'ztfLCDR1_%06d.hdf5';
DefV.ColCell            = {'HMJD','Mag','MagErr','ColorCoef','Flags'};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


FileName = sprintf('%s',InPar.FileTemp);
FileName = sprintf(FileName,FieldID);

Ncol = numel(InPar.ColCell);

Cat = HDF5.load(FileName,'/AllLC',[Lines(1) 1],[Lines(2)-Lines(1)+1, Ncol]);
ColCell = InPar.ColCell;

if (nargout>2)
    CatProp = HDF5.load(FileName,'/IndAllLC');
end