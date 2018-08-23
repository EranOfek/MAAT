function Cat=get_usnob1(RA,Dec,FOV,Shape);
%---------------------------------------------------------------------------
% get_usnob1 function                                             Catalogue
% Description: Get USNO-B1.0 catalog from local disk (using scat).
% Input  : - RA (J2000.0),    [HH MM SS] or [Radians]
%          - Dec (J2000.0), [Sign, DD MM SS] or [Radians]
%          - Half the Field of view (or radius), [arcsec]
%          - Field of view shape:
%            'b' - box, default.
%            'c' - circle.
% Output : - USNO-B1.0 catalog
%            [RA, Dec, MagB1  MagR1  MagB2  MagR2  MagN  PM  NI]
%            RA, Dec in radians.
% Comments: wcstools/scat + USNO-A2.0 should be installed
% Tested : MATLAB 5.3
%     By : Eran O. Ofek                October 2003
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
RAD = 180./pi;

if (nargin==3),
   Shape = 'b';
elseif (nargin==4),
   % do nothing
else
  error('Illegal number of input arguments');
end

if (length(RA)==1),
   RA = convertdms(RA,'r','H');
end
if (length(Dec)==1),
   Dec = convertdms(Dec,'r','D');
end

StrRA  = sprintf('%02d:%02d:%05.2f',RA);
if (Dec(1)>0),
   StrDec = sprintf('+%02d:%02d:%04.1f',Dec(2:4));
else
   StrDec = sprintf('-%02d:%02d:%04.1f',Dec(2:4));
end

TmpFile = 'tmp.matlab.usnob1';
eval(['!rm -rf ',TmpFile]);
switch Shape
 case 'b'
    RunStr1 = ['!setenv UB1_PATH /home/castor3/usnob; scat -c ub1 -d minpmqual=11 -s m -n 1000000 -r -',sprintf('%d',FOV),' -j ',StrRA,' ',StrDec,' >> ',TmpFile];
 case 'c'
    RunStr1 = ['!setenv UB1_PATH /home/castor3/usnob; scat -c ub1 -d minpmqual=11 -s m -n 1000000 -r ',sprintf('%d',FOV),' -j ',StrRA,' ',StrDec,' >> ',TmpFile];
 otherwise
    error('Unknown Shape option');
end

eval(RunStr1);

FID  = fopen(TmpFile,'r');
Line = fscanf(FID,'%f %f %f %f %f %f %f %f %d %d %f');
List = vec2mat(Line,11);
fclose(FID);

List(:,2:3) = List(:,2:3)./RAD;
Cat    = [List(:,[2, 3, 4, 5, 6, 7, 8, 9, 10])];
I      = find(Cat==99.99);
Cat(I) = NaN;

