function Cat=get_usnoa2(RA,Dec,FOV,OutputType,PathUSNO);
%---------------------------------------------------------------------------
% get_usnoa2 function                                             Catalogue
% Description: Get USNO-A2.0 catalog from local disk.
% Input  : - RA (J2000.0),    [HH MM SS] or [Radians], or sexagesimal string.
%          - Dec (J2000.0), [Sign, DD MM SS] or [Radians], or sexagesimal string.
%          - Half the Field of view, [arcsec]
%          - Output Type: 'dms' | 'rad'
%          - USNO-A2.0 path, default is '/home/castor3/USNO-A2.0'.
% Output : - USNO-A2.0 catalog
%            [RA_HH, RA_MM, RA_SS, Dec_Sign, Dec_DD, Dec_MM, Dec_SS, MagO, MagE]%            [RA_Rad, Dec_Rad, MagO, MagE]
%            for Output Type 'dms' and 'rad', respectively.
% Tested : MATLAB 5.3
%     By : Eran O. Ofek                October 2003
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Platform: Linux/UNIX only
% Needed : WCSTools scat program should be in the directory '../WCSTools/',
%          relative to the location of this program.
%          + USNO-A2.0 should be installed
%---------------------------------------------------------------------------
DirSCAT  = '../WCSTools';
Dir      = which_dir('get_usnoa2');
SCAT     = 'scat';
PathSCAT = sprintf('%s/%s/%s',Dir,DirSCAT,SCAT);


if (nargin==4),
   PathUSNO = '/home/castor3/USNO-A2.0';
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

RA  = convertdms(RA,'gH','H');
Dec = convertdms(Dec,'gD','D');

StrRA  = sprintf('%02d:%02d:%05.2f',RA);
if (Dec(1)>0),
   StrDec = sprintf('+%02d:%02d:%04.1f',Dec(2:4));
else
   StrDec = sprintf('-%02d:%02d:%04.1f',Dec(2:4));
end

TmpFile = 'tmp.matlab.usnoa2';
eval(['!rm -rf ',TmpFile]);
RunStr1 = ['!setenv UA2_PATH ',PathUSNO,'; ',PathSCAT,' -c ua2 -s m -n 10000 -r -',sprintf('%d',FOV),' -j ',StrRA,' ',StrDec,' >> ',TmpFile];


eval(RunStr1);

FID = fopen(TmpFile,'r');

I = 0;
while (feof(FID)==0),
   I         = I + 1;
   Line = fgetl(FID);
   if (Line(28)=='+'),
      Sign = 1;
   else
      Sign = -1;
   end
   RA   = [str2num(Line(15:16)), str2num(Line(18:19)), str2num(Line(21:26))];
   Dec  = [Sign, str2num(Line(29:30)), str2num(Line(32:33)), str2num(Line(35:39))];
   MagO = [str2num(Line(42:45))];
   MagE = [str2num(Line(48:51))];

   Cat(I,:)  = [RA, Dec, MagO, MagE];
end

switch OutputType
 case 'dms'
    % do nothing
 case 'rad'
    RA  = convertdms(Cat(:,1:3),'H','r');
    Dec = convertdms(Cat(:,4:7),'D','R');
    Cat = [RA, Dec, Cat(:,8:9)];
 otherwise
    error('Unknown OutputType option');
end

