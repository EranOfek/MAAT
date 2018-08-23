function wget_sdo_aia_image(DateStart)


DateStart = [2010 10 01]

JDstart = julday(DateStart([3 2 1]));

Itot = 0;
for I=1:1:700,
   JD = JDstart + I;
   Date = jd2date(JD);
   Date = Date([3 2 1]);


   Filter = 211;




   % Example: http://jsoc.stanford.edu/data/aia/synoptic/2010/05/13/H0000/AIA20100513_0000_0094.fits

   Type = 'jp2';
   
   switch lower(Type)
    case 'fits'
       URL = sprintf('http://jsoc.stanford.edu/data/aia/synoptic/%04d/%02d/%02d/H%02d00/',...
                     Date(1),Date(2),Date(3),Date(4));
       T=urlread(URL);
       MatchStr = sprintf('>AIA\\d\\d\\d\\d\\d\\d\\d\\d_\\d\\d\\d\\d_%04d.fits<',Filter);
       RE=regexp(T,MatchStr,'match');
    case 'jp2'
       URL = sprintf('http://jsoc.stanford.edu/data/aia/images/%04d/%02d/%02d/%d/',Date(1),Date(2),Date(3),Filter);
   
       T=urlread(URL);
       MatchStr = sprintf('>\\d\\d\\d\\d_\\d\\d_\\d\\d__\\d\\d_\\d\\d_\\d\\d_\\d\\d\\d__SDO_AIA_AIA_%d.jp2<',Filter);
       RE=regexp(T,MatchStr,'match');
    otherwise
       error('Unknown Type option');
   end
   
   
   Nim  = length(RE);
   Link = cell(Nim,1);
   for Iim=1:1:Nim,
      Image      = RE{Iim}(2:end-1);
      Link{Iim}  = sprintf('%s%s',URL,Image);
   end
   
   www.pwget(Link,'',15);
   
   
   
   Files = dir('*.jp2');
   Nf = length(Files);
   
   
   
   [MatX,MatY] = meshgrid([1:1:4096]',[1:1:4096]');
   MatR = sqrt((MatX - 2048).^2 + (MatY - 2048).^2);
   Iin = find(MatR<=1524);
   Ian = find(MatR>1524 & MatR<2048);
   
   Nbin = 4;
   Ncyc = floor(Nf./Nbin);
   Ncyc = Nf;
   for Icyc=1:1:Ncyc,
      Itot = Itot + 1;

      Iim = Icyc;
      %Im = zeros(4096,4096,3);
      %JDbin = zeros(Nbin,1);
      %for Ibin=1:1:Nbin,
      %   Iim = (Icyc-1).*Nbin+Ibin;
      %   Im(:,:,Ibin) = imread(Files(Iim).name);
      %   Date = datevec(Files(Iim).name,'yyyy_mm_dd__HH_MM_SS');
      %   JDbin(Ibin) = julday(Date([3 2 1 4 5 6]));
      %end
      %% bin data
      %MedIm = median(Im,3);
   
      MedIm    = imread(Files(Iim).name);
      Date  = datevec(Files(Iim).name,'yyyy_mm_dd__HH_MM_SS');
      JDbin = julday(Date([3 2 1 4 5 6]));
   
      Res(Itot).JD = mean(JDbin);
      Res(Itot).RangeJD = range(JDbin);
      Res(Itot).SumCount = sum(MedIm(:));
   
      Res(Itot).SumCountIn = sum(MedIm(Iin));
      Res(Itot).SumCountAn = sum(MedIm(Ian));

      clear MedIm

   end

   delete('*.jp2');
   save Res.mat Res

end


