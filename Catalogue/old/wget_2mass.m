function [Cat,ColCell,ColS]=wget_2mass(RA,Dec,SearchRadius)
%------------------------------------------------------------------------------
% wget_2mass function                                                Catalogue
% Description: A web based search of the 2MASS all sky point source catalog.
%              Installation: 1. install cdsclient (instructions can be found
%              in: http://cdsarc.u-strasbg.fr/doc/cdsclient.html).
%              2. edit the first line of code in this program to point
%              to the location of the find2mass script.
%              This program will work only from Linux.
%              This program is replacing search2mass.m
% Input  : - R.A. in [H M S] format or in radians or in sexagesimal string.
%          - Dec. in [Sign D M S] format or in radians or in sexagesimal
%            string.
%          - Search radius in arcsec.
% Output : - 2MASS point source within the searched area. Return -1 if failed.
%            Columns description:
%            1  - RA [deg]
%            2  - Dec [deg]
%            3  - J magnitude
%            4  - J mag error
%            5  - H magnitude
%            6  - H mag error
%            7  - K magnitude
%            8  - K mag error
%            9  - JD
%            10 - J mag chi2   
%            11 - H mag chi2   
%            12 - K mag chi2   
%            13 - distance in arcsec from USNO counterpart
%            14 - PA in deg from USNO counterpart
%            15 - USNO counterpart B mag
%            16 - USNO counterpart R mag
%          - Cell array containing the table column names.
%          - Structure containing the table column names.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      July 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Needed : The exacutable program find2mass should be
%          in the same directory as search2mass.m program
% Example: Cat=wget_2mass([10 0 0],[+1 40 0 0],10);
% Reliable: 1
%------------------------------------------------------------------------------
SearchProgExe = 'bin/vizquery/cdsclient-3.71/find2mass';
Dir = which_dir('search2mass');
SEARCH_PROGRAM = sprintf('%s/../%s',Dir,SearchProgExe);

TMP_CAT        = tempname;
RAD  = 180./pi;  

RA  = convertdms(RA,'gH','r');
Dec = convertdms(Dec,'gD','R');

% maximum number of source is 100000
RunStr = sprintf('%s -c %f %f -r %f -max 200000 > %s',SEARCH_PROGRAM,RA.*RAD,Dec.*RAD,SearchRadius./60,TMP_CAT);

RunStr
system(RunStr);
system(sprintf('sed "s#|# #g" %s > %s_1',TMP_CAT,TMP_CAT));
system(sprintf('mv %s_1 %s',TMP_CAT,TMP_CAT));

%--- read 2MASS catalog ---
FID  = fopen(TMP_CAT,'r');
if (FID==-1),
   % file not exist
   Cat = -1;
else

   Line = fgetl(FID);
   Line = fgetl(FID);
   Line = fgetl(FID);
   Line = fgetl(FID);
   Line = fgetl(FID);
      
   NL = wc(TMP_CAT,'l')-6;
   %Cat  = zeros(NL,16);
   Cat  = zeros(NL,8);
   for LineI=1:1:NL,
      ObjRA     = str2num_nan(fscanf(FID,'%s',1));
      ObjDec    = str2num_nan(fscanf(FID,'%s',1));
      ErrMaj    = str2num_nan(fscanf(FID,'%s',1));
      ErrMin    = str2num_nan(fscanf(FID,'%s',1));
      ErrPA     = str2num_nan(fscanf(FID,'%s',1));
      Name      = str2num_nan(fscanf(FID,'%s',1));
      Jmag      = str2num_nan(fscanf(FID,'%s',1));
      JmagErr   = str2num_nan(fscanf(FID,'%s',1));
      JmagME    = str2num_nan(fscanf(FID,'%s',1));
      JmagSN    = str2num_nan(fscanf(FID,'%s',1));
      Hmag      = str2num_nan(fscanf(FID,'%s',1));
      HmagErr   = str2num_nan(fscanf(FID,'%s',1));
      HmagME    = str2num_nan(fscanf(FID,'%s',1));
      HmagSN    = str2num_nan(fscanf(FID,'%s',1));
      Kmag      = str2num_nan(fscanf(FID,'%s',1));
      KmagErr   = str2num_nan(fscanf(FID,'%s',1));
      KmagME    = str2num_nan(fscanf(FID,'%s',1));
      KmagSN    = str2num_nan(fscanf(FID,'%s',1));

      Line = fgetl(FID);

if (0==1),
      Qfl       = fscanf(FID,'%s',1);
      Rfl       = fscanf(FID,'%s',1);
      Bfl       = fscanf(FID,'%s',1);
      Cfl       = fscanf(FID,'%s',1);
      N_images  = fscanf(FID,'%s',1);
      Prox      = fscanf(FID,'%s',1);
      ProxPA    = fscanf(FID,'%s',1);
      ProxCen   = fscanf(FID,'%s',1);
      X         = fscanf(FID,'%s',1);
      A         = fscanf(FID,'%s',1);
      Cntr      = fscanf(FID,'%d',1);
      H         = fscanf(FID,'%s',1);
      ObsDate   = fscanf(FID,'%s',1);
      Scan      = fscanf(FID,'%s',1);
      Glon      = fscanf(FID,'%f',1);
      Glat      = fscanf(FID,'%f',1);
      Xscan     = fscanf(FID,'%f',1);
      JD        = fscanf(FID,'%f',1);
      Jchi      = str2num_nan(fscanf(FID,'%s',1));
      Hchi      = str2num_nan(fscanf(FID,'%s',1));
      Kchi      = str2num_nan(fscanf(FID,'%s',1));
      Japm      = fscanf(FID,'%s',1);
      JapmErr   = fscanf(FID,'%s',1);
      Hapm      = fscanf(FID,'%s',1);
      HapmErr   = fscanf(FID,'%s',1);
      Kapm      = fscanf(FID,'%s',1);
      KapmErr   = fscanf(FID,'%s',1);
      EdgeN     = fscanf(FID,'%s',1);
      EW        = fscanf(FID,'%s',1);
      QQQ       = fscanf(FID,'%s',1);
      D         = fscanf(FID,'%s',1);
      U         = fscanf(FID,'%s',1);
      O         = fscanf(FID,'%s',1);
      USNO_Rad  = str2num_nan(fscanf(FID,'%s',1));
      USNO_PA   = str2num_nan(fscanf(FID,'%s',1));
      USNO_Bmag = str2num_nan(fscanf(FID,'%s',1));
      USNO_Rmag = str2num_nan(fscanf(FID,'%s',1));
      N         = fscanf(FID,'%d',1);
      Extnd     = fscanf(FID,'%s',1);
      ScanN     = fscanf(FID,'%s',1);
      CoAdded   = fscanf(FID,'%s',1);
      Num       = fscanf(FID,'%s',1);
      Junk      = fscanf(FID,'%s',1);
end

      Cat(LineI,:) = [ObjRA, ObjDec, Jmag, JmagErr, Hmag, HmagErr, Kmag, KmagErr];
      %, JD, Jchi, Hchi, Kchi, USNO_Rad, USNO_PA, USNO_Bmag, USNO_Rmag];
   end
   
   fclose(FID);
   
   delete(TMP_CAT);
end


if (nargout>1),
   ColCell = {'RA','Dec','J','JErr','H','HErr','K','KErr','JD','Jchi','Hchi','Kchi','USNOrad','USNOpa','USNO_B','USNO_R'};
   ColS = cell2struct(num2cell([1:1:length(ColCell)]),ColCell,2)
end



function Number=test_str(String);
if (String(1)=='-'),
   Number = NaN;
else
   Number = str2num(String);
end

delete(TMP_CAT);
