function [Cat,ColCell,Col,UnitsCell]=wget_2mass(RA,Dec,SearchSize,SearchShape)
%------------------------------------------------------------------------------
% wget_2mass function                                                Catalogue
% Description: Query the 2MASS catalog using the VizieR web service.
% Installation: 1. install cdsclient (instructions can be found
%               in: http://cdsarc.u-strasbg.fr/doc/cdsclient.html)
%               in $USER/matlab/fun/bin/vizquery/cdsclient-3.71/
%               2. If you installed the cdsclient in a different location,
%               then edit the first few lines of the code accordingly.
%               This program is replacing search2mass.m
% Input  : - R.A. in [H M S] format or in radians or in sexagesimal string.
%          - Dec. in [Sign D M S] format or in radians or in sexagesimal
%            string.
%          - Search radius in radians.
%          - Search shape {'circle','box'}, default is 'circle'.
% Output : - Matrix of 2MASS sources found within search region.
%            The matrix contains the following columns:
%            1. RA (J2000) [rad].
%            2. Dec (J2000) [rad].
%            3. J mag [mag].
%            4. J ma err [mag].
%            5. H mag [mag].
%            6. H ma err [mag].
%            7. K mag [mag].
%            8. K ma err [mag].
%            9. JD of observation [days].
%            10. Distance to nearest USNO-A2 source.
%            11. PA of nearst USNO-A2 source.
%            12. B mag of nearest source [mag].
%            12. R mag of nearest source [mag].
%          - Cell array of column names.
%          - A structure in which the fields are the column names,
%            and the field value is the column index.
%          - Cell array of column units.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jul 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Needed : The exacutable program find2mass should be
%          in the same directory as search2mass.m program
% Example: [Cat,ColCell,Col,Units]=wget_2mass([10 0 0],[+1 40 0 0],100./(RAD.*3600));
% Reliable: 2
%------------------------------------------------------------------------------
RAD  = 180./pi;
MaxNumStar  = 100000;

ProgName = 'find2mass';
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


%sprintf('%s %s -c %s %s %7.1f -m %d',Prog, RA, Dec, ShapeKey, SearchSize, MaxNumStar);
SearchSizeAS = SearchSize.*RAD.*3600;  % search size in arcsec
[Stat,Res]=system(sprintf('%s -c %s %s %s %7.1f -m %d',Prog, RA, Dec, ShapeKey, SearchSizeAS, MaxNumStar));

[Lines]=textscan(Res,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n','commentstyle','#','delimiter','|');

Nl = length(Lines{1});
if (Nl==0),
  Out = [];
end

ColCell    = {'RA','Dec','J','errJ','H','errH','K','errK','JD','Dist','PA','B','R'};
Cat        = zeros(Nl,length(ColCell)).*NaN;
Col        = cell2struct(num2cell([1:1:length(ColCell)]),ColCell,2);
UnitsCell  = {'rad','rad','mag','mag','mag','mag','mag','mag','day','arcsec','deg','mag','mag'};


Nl = length(Lines{1});
for Il=1:1:Nl,
   Info                  = textscan(Lines{1}{Il},'%f %f');
   Cat(Il,1) = Info{1}./RAD;
   Cat(Il,2) = Info{2}./RAD;

   Info                  = textscan(Lines{4}{Il},'%f %f %f %f');
   Cat(Il,3) = str2num_nan(Info{1});
   Cat(Il,4) = str2num_nan(Info{2});

   Info                  = textscan(Lines{5}{Il},'%f %f %f %f');
   Cat(Il,5) = str2num_nan(Info{1});
   Cat(Il,6) = str2num_nan(Info{2});

   Info                  = textscan(Lines{6}{Il},'%f %f %f %f');
   Cat(Il,7) = str2num_nan(Info{1});
   Cat(Il,8) = str2num_nan(Info{2});

   Info                  = textscan(Lines{12}{Il},'%f %f');
   Cat(Il,9) = Info{2};  % ObsJD 

   Info                  = textscan(Lines{19}{Il},'%s %f %f');
   Cat(Il,10) = str2num_nan(Info{2});
   Cat(Il,11) = str2num_nan(Info{3});
   
   Info                  = textscan(Lines{20}{Il},'%f %f %f %*[^\n]');
   Cat(Il,12) = str2num_nan(Info{1});
   Cat(Il,13) = str2num_nan(Info{2});

end

% sort by declination
Cat = sortrows(Cat,Col.Dec);
