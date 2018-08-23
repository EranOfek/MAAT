function [Cat,ColCell,Col,UnitsCell]=wget_usnob1(RA,Dec,SearchSize,SearchShape)
%------------------------------------------------------------------------------
% wget_usnob1 function                                               Catalogue
% Description: Query the USNO-B1 catalog using the VizieR web service.
% Installation: 1. install cdsclient (instructions can be found
%               in: http://cdsarc.u-strasbg.fr/doc/cdsclient.html)     
%               in $USER/matlab/fun/bin/vizquery/cdsclient-3.4/
%               2. If you installed the cdsclient in a different location,
%               then edit the first few lines of the code accordingly.
%               This program is replacing search2mass.m
% Input  : - J2000.0 R.A. in radians, [H M S] or sexagesimal string.
%          - J2000.0 Dec. in radians, [Sign D M S] or sexagesimal string.
%          - Search radius or search box size in radians.
%          - Search shape {'circle','box'}, default is 'circle'.
% Output : - Matrix of USNO-B1 sources found within search region.
%            The matrix contains the following columns:
%            1. RA (J2000) [rad].
%            2. Dec (J2000) [rad].
%            3. err in RA [mas].
%            4. err in Dec [mas].
%            5. PM in RA [mas/yr].
%            6. PM in Dec [mas/yr].
%            7. Number of detections.
%            8. B1mag [mag].
%            9. R1mag [mag].
%            10. B2mag [mag].
%            11. R2mag [mag].
%            12. I1mag [mag].
%          - Cell array of column names.
%          - A structure in which the fields are the column names,
%            and the field value is the column index.
%          - Cell array of column units.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jun 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cat]=wget_usnob1(1,1,30./(3600.*RAD));
%          [Cat,ColCell,Col,UnitsCell] = wget_usnob1('14:56:01.1',[-1 67 01 10],30./(3600.*RAD));
% Reliable: 2
%------------------------------------------------------------------------------
RAD  = 180./pi;
MaxNumStar  = 100000;

ProgName = 'findusnob1';
Path     = vizquery_path(ProgName);   % GO THIS THIS PROGRAM AND EDIT IF NEEDED
Prog     = sprintf('%s%s',Path,ProgName);

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

%sprintf('%s %s %s %s %7.1f -m %d',Prog, RA, Dec, ShapeKey, SearchSize, MaxNumStar)
SearchSizeAS = SearchSize.*RAD.*3600;  % search size in arcsec
[Stat,Res]=system(sprintf('%s %s %s %s %7.1f -m %d',Prog, RA, Dec, ShapeKey, SearchSizeAS, MaxNumStar));


[Lines]=textscan(Res,'%s %s %s %s %s %s %s\n','commentstyle','#','delimiter','|');


Nl = length(Lines{1});
if (Nl==0),
   Cat = [];
end

ColCell    = {'RA','Dec','sRA','sDec','pmRA','pmDec','Ndet','B1mag','R1mag','B2mag','R2mag','Imag'};
Cat        = zeros(Nl,length(ColCell)).*NaN;
Col        = cell2struct(num2cell((1:1:length(ColCell))),ColCell,2);
UnitsCell  = {'rad','rad','mas','mas','mas/yr','mas/yr','','mag','mag','mag','mag','mag'};


for Il=1:1:Nl,
   Cat(Il,Col.RA)        = convertdms(Lines{1}{Il}(1:12),'SH','r');
   Cat(Il,Col.Dec)       = convertdms(Lines{1}{Il}(13:24),'SD','R');
   Info                  = textscan(Lines{1}{Il},'%s %f %f %f %f %s %f %f %f %f %f %s');
   Cat(Il,Col.sRA)       = Info{2};   % err in RA [mas]
   Cat(Il,Col.sDec)      = Info{3};   % err in Dec [mas]
   Cat(Il,Col.pmRA)      = Info{4};   % PM in RA [mas/yr]
   Cat(Il,Col.pmDec)     = Info{5};   % PM in Dec [mas/yr]
   Cat(Il,Col.Ndet)      = Info{11};

   Cat(Il,Col.B1mag) = str2num_nan(Info{1}{1});
   Info              = textscan(Lines{3}{Il},'%s');
   Cat(Il,Col.R1mag) = str2num_nan(Info{1}{1});
   Info              = textscan(Lines{4}{Il},'%s');
   Cat(Il,Col.B2mag) = str2num_nan(Info{1}{1});
   Info              = textscan(Lines{5}{Il},'%s');
   Cat(Il,Col.R2mag) = str2num_nan(Info{1}{1});
   Info              = textscan(Lines{6}{Il},'%s');
   Cat(Il,Col.Imag)  = str2num_nan(Info{1}{1});

end   

% sort by declination
Cat = sortrows(Cat,Col.Dec);
