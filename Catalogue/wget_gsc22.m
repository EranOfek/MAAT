function Catalog=wget_gsc22(RA,Dec,Radius,MagLimit);
%-----------------------------------------------------------------------------
% wget_gsc22 function                                               Catalogue
% Description: Get GSC-2.2 catalog from STSCI web server
% Input  : - RA, [H M S] or [Radians], or sexagesimal string.
%          - Dec [sign D M S] or [Radians], or sexagesimal string.
%          - Search radius in arcmin, default is 12.
%          - Magnitude limit, default is 21.
% Output : - Catalog with the following columns:
%            RA    = 1
%            Dec   = 2
%            RAe   = 3
%            Dece  = 4
%            Epoch = 5
%            RAPM  = 6
%            DecPM = 7
%            RAPMe = 8
%            DecPMe= 9
%            Fmag  = 10
%            FmagE = 11
%            Jmag  = 12
%            JmagE = 13
%            Vmag  = 14
%            VmagE = 15
%            Nmag  = 16
%            NmagE = 17
%            A     = 18
%            E     = 19
%            PA    = 20
%            C     = 21
%            Status= 22
%            N     = 22
% Reference: http://www-gsss.stsci.edu/support/data_access.htm
% Tested : Matlab 5.3
%     By : Eran O. Ofek            June 2003
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: C=wget_gsc22([10 0 0],[1 40 0 0],10,25);
%-----------------------------------------------------------------------------
if (nargin==2),
   Radius = 12;
   MagLimit = 21;
elseif (nargin==3),
   MagLimit = 21;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end


RA  = convertdms(RA,'gH','H');
Dec = convertdms(Dec,'gD','D');


MaxN     = 10000;

ColRA    = 1;
ColDec   = 2;
ColRAe   = 3;
ColDece  = 4;
ColEpoch = 5;
ColRAPM  = 6;
ColDecPM = 7;
ColRAPMe = 8;
ColDecPMe= 9;
ColFmag  = 10;
ColFmagE = 11;
ColJmag  = 12;
ColJmagE = 13;
ColVmag  = 14;
ColVmagE = 15;
ColNmag  = 16;
ColNmagE = 17;
ColA     = 18;
ColE     = 19;
ColPA    = 20;
ColC     = 21;
ColStatus= 22;
ColN     = 22;

!rm gsc22query*
!rm tmp.Catalog

if (Dec(1)==-1),
   StrURL = sprintf('"http://www-gsss.stsci.edu/cgi-bin/gsc22query.exe?ra=%02d+%02d+%02d&dec=-%02d+%02d+%02d&r1=0&r2=%07.3f&m1=0.0&m2=%05.1f&n=%d&submit2=Submit+Request"',RA(1), RA(2), RA(3), Dec(2), Dec(3), Dec(4), Radius, MagLimit, MaxN)
elseif (Dec(1)==1),
   StrURL = sprintf('"http://www-gsss.stsci.edu/cgi-bin/gsc22query.exe?ra=%02d+%02d+%02d&dec=+%02d+%02d+%02d&r1=0&r2=%07.3f&m1=0.0&m2=%05.1f&n=%d&submit2=Submit+Request"',RA(1), RA(2), RA(3), Dec(2), Dec(3), Dec(4), Radius, MagLimit, MaxN)
else
   error('Unkown Dec format');
end
eval (['!wget ',StrURL]);

%Rm = ls('gsc22query*');
!mv gsc22query* tmp.Catalog

FID = fopen('tmp.Catalog','r');
Line     = fgetl(FID);
Line     = fgetl(FID);

Catalog = zeros(0,ColN);

I       = 0;
Name    = fscanf(FID,'%s',1);
while (Name(1)=='S' | Name(1)=='N');
   I = I + 1;
   Catalog(I,1)  = fscanf(FID,'%f',1);
   Catalog(I,2)  = fscanf(FID,'%f',1);
   Catalog(I,3)  = fscanf(FID,'%f',1);
   Catalog(I,4)  = fscanf(FID,'%f',1);
   Catalog(I,5)  = fscanf(FID,'%f',1);
   Catalog(I,6)  = fscanf(FID,'%f',1);
   Catalog(I,7)  = fscanf(FID,'%f',1);
   Catalog(I,8)  = fscanf(FID,'%f',1);
   Catalog(I,9)  = fscanf(FID,'%f',1);
   Catalog(I,10) = fscanf(FID,'%f',1);
   Catalog(I,11) = fscanf(FID,'%f',1);
   Catalog(I,12) = fscanf(FID,'%f',1);
   Catalog(I,13) = fscanf(FID,'%f',1);
   Catalog(I,14) = fscanf(FID,'%f',1);
   Catalog(I,15) = fscanf(FID,'%f',1);
   Catalog(I,16) = fscanf(FID,'%f',1);
   Catalog(I,17) = fscanf(FID,'%f',1);
   Catalog(I,18) = fscanf(FID,'%f',1);
   Catalog(I,19) = fscanf(FID,'%f',1);
   Catalog(I,20) = fscanf(FID,'%f',1);
   Catalog(I,21) = fscanf(FID,'%f',1);
   Catalog(I,22) = fscanf(FID,'%f',1);
   Name          = fscanf(FID,'%s',1);
end

fclose(FID);

% remove catalog
!rm tmp.Catalog
