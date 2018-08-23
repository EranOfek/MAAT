function [Cat,Header]=search2mass(RA,Dec,SearchRadius);
%----------------------------------------------------------------------------
% search2mass function                                             Catalogue
% Description: A web based search of the 2MASS all sky point source catalog.
%              The search is done using the find2mass utility and aclient.
% Input  : - R.A. in [H M S] format or in radians or in sexagesimal string.
%          - Dec. in [Sign D M S] format or in radians or in sexagesimal string.
%          - Search radius in arcmin.
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
%          - Cell array containing the table columns.
% Tested : Matlab 5.3
%     By : Eran O. Ofek              July 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Needed : The exacutable program find2mass should be
%          in the same directory as search2mass.m program
% Example: Cat=search2mass([10 0 0],[+1 40 0 0],10);
%----------------------------------------------------------------------------
Dir = which_dir('search2mass');
SearchProgExe = 'find2mass';
SEARCH_PROGRAM = sprintf('%s/%s',Dir,SearchProgExe);

TMP_CAT        = 'tmp_2mass.txt';
RAD  = 180./pi;  

RA  = convertdms(RA,'gH','r');
Dec = convertdms(Dec,'gD','R');

% maximum number of source is 100000
RunStr = sprintf('%s %f %f -r %f -m 100000 > %s',SEARCH_PROGRAM,RA.*RAD,Dec.*RAD,SearchRadius,TMP_CAT);

RunStr
system(RunStr);


%--- read 2MASS catalog ---
FID  = fopen(TMP_CAT,'r');
if (FID==-1),
   % file not exist
   Cat = -1;
else

   Line = fgetl(FID);
   Line = fgetl(FID);
   Line = fgetl(FID);
   
   
   NL = wc(TMP_CAT,'l')-4;
   Cat  = zeros(NL,16);
   for LineI=1:1:NL,
      ObjRA     = fscanf(FID,'%f',1);
      ObjDec    = fscanf(FID,'%f',1);
      ErrMaj    = fscanf(FID,'%f',1);
      ErrMin    = fscanf(FID,'%f',1);
      ErrPA     = fscanf(FID,'%f',1);
      Name      = fscanf(FID,'%s',1);
      Jmag      = test_str(fscanf(FID,'%s',1));
      JmagErr   = test_str(fscanf(FID,'%s',1));
      JmagME    = test_str(fscanf(FID,'%s',1));
      JmagSN    = test_str(fscanf(FID,'%s',1));
      Hmag      = test_str(fscanf(FID,'%s',1));
      HmagErr   = test_str(fscanf(FID,'%s',1));
      HmagME    = test_str(fscanf(FID,'%s',1));
      HmagSN    = test_str(fscanf(FID,'%s',1));
      Kmag      = test_str(fscanf(FID,'%s',1));
      KmagErr   = test_str(fscanf(FID,'%s',1));
      KmagME    = test_str(fscanf(FID,'%s',1));
      KmagSN    = test_str(fscanf(FID,'%s',1));
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
      Jchi      = test_str(fscanf(FID,'%s',1));
      Hchi      = test_str(fscanf(FID,'%s',1));
      Kchi      = test_str(fscanf(FID,'%s',1));
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
      USNO_Rad  = test_str(fscanf(FID,'%s',1));
      USNO_PA   = test_str(fscanf(FID,'%s',1));
      USNO_Bmag = test_str(fscanf(FID,'%s',1));
      USNO_Rmag = test_str(fscanf(FID,'%s',1));
      N         = fscanf(FID,'%d',1);
      Extnd     = fscanf(FID,'%s',1);
      ScanN     = fscanf(FID,'%s',1);
      CoAdded   = fscanf(FID,'%s',1);
      Num       = fscanf(FID,'%s',1);
      Junk      = fscanf(FID,'%s',1);

      Cat(LineI,:) = [ObjRA, ObjDec, Jmag, JmagErr, Hmag, HmagErr, Kmag, KmagErr, JD, Jchi, Hchi, Kchi, USNO_Rad, USNO_PA, USNO_Bmag, USNO_Rmag];
   end
   
   fclose(FID);
   
   delete(TMP_CAT);
end


if (nargout>1),
   Header = {'RA','Dec','J','JErr','H','HErr','K','KErr','JD','Jchi','Hchi','Kchi','USNOrad','USNOpa','USNO_B','USNO_R'};
end



function Number=test_str(String);
if (String(1)=='-'),
   Number = NaN;
else
   Number = str2num(String);
end
