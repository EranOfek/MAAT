function MergedData=gaia_dr1_readall_select(DecRange)
% Select GAIA sources in Dec zone for constructing GAIA catalog
% Pacakge: VO.prep.GAIA


DataDir = '/raid/eran/Catalogue/GAIA-DR1/Orig';

PWD = pwd;

cd(DataDir);
Files = dir('*.hdf5');
Nf    = numel(Files);

MergedData = zeros(0,6);
for If=1:1:Nf
    sprintf('File %d out of %d\r',If,Nf);
    Data = Util.IO.loadh(Files(If).name);
    Flag = Data(:,3)>=DecRange(1) & Data(:,3)<=DecRange(2);
    MergedData = [MergedData; Data(Flag,:)];
end
