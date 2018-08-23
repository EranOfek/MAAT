function Count=gaia_dr1_readall2hdf5
% Create an HDF5 version of the GAIA-DR1 files with a subset of columns.
% Example: Count=VO.prep.GAIA.gaia_dr1_readall2hdf5

DataDir = '/raid/eran/Catalogue/GAIA-DR1/Orig';

PWD = pwd;

cd(DataDir);
Files = dir('*.csv');
Nf    = numel(Files);

Count = 0;
for If=1:1:Nf
    sprintf('File %d out of %d\r',If,Nf);
    [Data,ColCell,C,Flag]=VO.prep.GAIA.gaia_dr1_read_file(Files(If).name);
    FileHDF5 = sprintf('%s.hdf5',Files(If).name(1:end-4));
    
    %saveh(FileHDF5,Data(Flag,:));   % only good data
    saveh(FileHDF5,Data(:,:));   % keep all data
    Count = Count + sum(Flag);
end

% count:  664571170 