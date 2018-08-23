function [Out]=wget_usnob1(RA,Dec,SearchSize,SearchShape);
%------------------------------------------------------------------------------
% wget_usnob1 function                                               Catalogue
% Description: Query the USNO-B1 catalog using the VizieR web service.
%              Instellation: 
% Input  : - J2000.0 R.A. in radians, [H M S] or sexagesimal string.
%          - J2000.0 Dec. in radians, [Sign D M S] or sexagesimal string.
%          - Search radius or search box size in arcsec.
%          - Search shape {'circle','box'}, default is 'circle'.
% Output : - Structure containing the USNO-B1 catalog.
%            The structure contains the following fields:
%            .Name
%            .RA
%            .Dec
%            .sRA
%            .sDec
%            .Epoch
%            .pmRA
%            .pmDec
%            .Ndet
%            .B1mag
%            .B1sg
%            .R1mag
%            .R1sg
%            .B2mag
%            .B2sg
%            .R2mag
%            .R2sg
%            .I1mag
%            .I1sg
%            .Dist
% Tested : Matlab 7.3
%     By : Eran O. Ofek                         June 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [Out]=wget_usnob1(1,1,30);
%------------------------------------------------------------------------------
RAD  = 180./pi;
MaxNumStar  = 100000;

Dir = which_dir('wget_usnob1');
ProgPath = '../bin/vizquery/cdsclient-3.4/';
ProgName = 'findusnob1';
Prog     = sprintf('%s/%s%s',Dir,ProgPath,ProgName);


DefSearchShape   = 'circle';
if (nargin==3),
   SearchShape   = DefSearchShape;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end

RA  = convertdms(RA,'gH','SH');
Dec = convertdms(Dec,'gD','SD');

switch SearchShape
 case 'circle'
    ShapeKey = '-rs';
 case 'box'
    ShapeKey = '-bs';
 otherwise
    error('Unknown SearchShape option');
end

sprintf('%s %s %s %s %7.1f -m %d',Prog, RA, Dec, ShapeKey, SearchSize, MaxNumStar)
[Stat,Res]=system(sprintf('%s %s %s %s %7.1f -m %d',Prog, RA, Dec, ShapeKey, SearchSize, MaxNumStar));


[Lines]=strread(Res,'%263c%*[^\n]','commentstyle','shell');

Nl = size(Lines,1);
if (Nl==0),
   Out = [];
end
size(Lines)
for Il=1:1:Nl,
   Out.Name{Il,1}    = Lines(Il,1:12);
   Out.RA(Il,1)      = str2num_nan(Lines(Il,27:36))./RAD;
   Out.Dec(Il,1)     = str2num_nan(Lines(Il,37:46))./RAD;
   Out.sRA(Il,1)     = str2num_nan(Lines(Il,48:50));
   Out.sDec(Il,1)    = str2num_nan(Lines(Il,52:54));
   Out.Epoch(Il,1)   = str2num_nan(Lines(Il,56:61));
   Out.pmRA(Il,1)    = str2num_nan(Lines(Il,63:68));
   Out.pmDec(Il,1)   = str2num_nan(Lines(Il,70:75));
   Out.Ndet(Il,1)    = str2num_nan(Lines(Il,91));
   Out.B1mag(Il,1)   = str2num_nan(Lines(Il,98:102));
   Out.B1sg(Il,1)    = str2num_nan(Lines(Il,112:113));   % 0-non-stellar; 11-stellar
   Out.R1mag(Il,1)   = str2num_nan(Lines(Il,129:133));
   Out.R1sg(Il,1)    = str2num_nan(Lines(Il,143:144));  
   Out.B2mag(Il,1)   = str2num_nan(Lines(Il,160:164));
   Out.B2sg(Il,1)    = str2num_nan(Lines(Il,174:175));  
   Out.R2mag(Il,1)   = str2num_nan(Lines(Il,191:195));
   Out.R2sg(Il,1)    = str2num_nan(Lines(Il,205:206));  
   Out.I1mag(Il,1)   = str2num_nan(Lines(Il,222:226));
   Out.I1sg(Il,1)    = str2num_nan(Lines(Il,236:237));  
   Out.Dist(Il,1)    = str2num_nan(Lines(Il,255:263)); 
end   
