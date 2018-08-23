function Res=ccf_with_gaia(Spec);

Dir = which_dir('ccf_with_gaia');
GaiaDir = '../../data/AstroSpec/GAIA/Spec1A';

TempDir = dir(sprintf('%s%sT*',GaiaDir,filesep));
Ntd     = length(TempDir);

SmoothGaia = 5;
RedshiftVec = [-0.01:0.0001:0.01]';


GaiaWave = [2501:1:10500]';


Res = [];
Counter = 0;
for Itd=1:1:Ntd,
   CurDir   = sprintf('%s%s%s%s',GaiaDir,filesep,TempDir(Itd).name,filesep);
   FileDir  = dir(sprintf('%s*V000K1*.ASC',CurDir));
   Nf       = length(FileDir);
   for If=1:1:Nf,
      Counter = Counter + 1;

      File = sprintf('%s%s',CurDir,FileDir(If).name);

      GaiaSpecFlux = load(File);
      GaiaSpec = [GaiaWave, GaiaSpecFlux];

      GaiaSpec(:,2) = medfilt1(GaiaSpec(:,2),SmoothGaia);

      [CC,Npt,RedshiftVec]=ccf_spectra(Spec,GaiaSpec,RedshiftVec,'none','median');
      
      %[Alpha,RMS,Best]=chi2_spectra(Spec,GaiaSpec,RedshiftVec);
      %[MinRMS,MinInd] = min(RMS);
      [MaxCC,MaxInd] = max(CC);

      Res(Counter).MaxCC   = MaxCC;
      Res(Counter).BestZ    = RedshiftVec(MaxInd);

      Res(Counter).Name = FileDir(If).name;
      Res(Counter).FullPath = File;
      %Res(Counter).MinRMS   = MinRMS;
      %Res(Counter).BestZ    = RedshiftVec(MinInd);

      Res(Counter).Temp   = str2num(FileDir(If).name(2:6));
      Res(Counter).Grav   = str2num(FileDir(If).name(8:9));
      if (strcmp(FileDir(If).name(10),'M')==1),
         MetaSign = -1;
       else
         MetaSign = 1;
      end
      Res(Counter).Meta   = MetaSign.*str2num(FileDir(If).name(11:12));
      Res(Counter).Velo   = str2num(FileDir(If).name(14:16));

Counter
   end
end

