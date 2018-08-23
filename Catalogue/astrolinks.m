function [URL,Link]=astrolinks(Data,ColCoo,OutputFile,varargin);
%------------------------------------------------------------------------------
% astrolinks function                                                Catalogue
% Description: Given a list of coordinates, return the URL links to the
%              following webpages: SDSS chart, SDSS object, NED, FIRST, NVSS, 
%              and DSS. Moreover, generate an html page with the catalog data 
%              and links to the verious web pages.
% Input  : - Matrix or cell array. The matrix should contain
%            J2000.0 RA and Dec coordinates in decimal degrees.
%          - Indices of columns in matrix/cell array that contains
%            the coordiantes [RA, Dec].
%          - Name of output html file, default: no output file.
%          * arbitrary number of pairs of input arguments
%            (...,keyword,value,...) to be passed to html_table.m
%            See html_table.m for help.
% Output : - Cell array of links structure for each target.
%            URL.SDSSc, URL.SDSSo
%            URL.NED
%            URL.FIRST
%            URL.NVSS
%            URL.SIMBAD
%            URL.DSS2R, URL.DSS2B, URL.DSS2I, URL.DSS1R, URL.DSS1B
%          - Cell array of strings of links
% Tested : Matlab 7.0
%     By : Eran O. Ofek                   October 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;

NED_Radius    = 2;     % NED search radius in arcmin
DSS_FOV       = 12;    % DSS FOV in arcmin
SIMBAD_Radius = 0.5;

if (nargin==2),
   OutputFile = [];
else
   % do nothing
end



if (iscell(Data)==0),
   DataCell = num2cell(Data);
else
   DataCell = Data;
end
N      = size(DataCell,1);
RA     = zeros(N,1);
Dec    = zeros(N,1);
StrRA  = cell(N,1);
StrDec = cell(N,1);

% convert coordiantes and set links
for I=1:1:N,
   RA(I)     = DataCell{I,ColCoo(1)};            % deg
   Dec(I)    = DataCell{I,ColCoo(2)};            % deg
   StrRA{I}  = convertdms(RA(I),'d','SH');       % sexagesimal
   StrDec{I} = convertdms(Dec(I),'d','SD');      % sexagesimal
   
   % link to SDSS chart
   Link.SDSSc{I} = sprintf('<a href="http://cas.sdss.org/astro/en/tools/chart/navi.asp?ra=%10.6f&dec=%10.6f&scale=0.2&width=512&height=512&opt=&query=" target="_blank">SDSSc</a>',RA(I),Dec(I));
   URL.SDSSc{I} = sprintf('http://cas.sdss.org/astro/en/tools/chart/navi.asp?ra=%10.6f&dec=%10.6f&scale=0.2&width=512&height=512&opt=&query=',RA(I),Dec(I));

   % link to SDSS object
   Link.SDSSo{I} = sprintf('<a href="http://cas.sdss.org/astro/en/tools/explore/obj.asp?ra=%10.6f&dec=%10.6f" target="_blank">SDSSo</a>',RA(I),Dec(I));
   URL.SDSSo{I} = sprintf('http://cas.sdss.org/astro/en/tools/explore/obj.asp?ra=%10.6f&dec=%10.6f',RA(I),Dec(I));

   % Link to NED
   Link.NED{I} = sprintf('<a href="http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon=%10.6fd&lat=%10.6fd&radius=%05.2f&search_type=Near+Position+Search&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=500&img_stamp=YES&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY" target="_blank">NED</a>',RA(I),Dec(I),NED_Radius);
   URL.NED{I} = sprintf('http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon=%10.6fd&lat=%10.6fd&radius=%05.2f&search_type=Near+Position+Search&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=500&img_stamp=YES&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY',RA(I),Dec(I),NED_Radius);

   % Link to SIMBAD
   Link.SIMBAD{I} = sprintf('<a href="http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=J2000&CooEpoch=2000&Coord=%11.6f%%20%11.6f&submit=submit%%20query&Radius.unit=arcmin&CooEqui=2000&CooFrame=ICRS&Radius=%f">SIMBAD</a>',RA(I),Dec(I),SIMBAD_Radius);
   URL.SIMBAD{I} = sprintf('http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=J2000&CooEpoch=2000&Coord=%11.6f%%20%11.6f&submit=submit%%20query&Radius.unit=arcmin&CooEqui=2000&CooFrame=ICRS&Radius=%f',RA(I),Dec(I),SIMBAD_Radius);


   % Link to FIRST
   Link.FIRST{I} = sprintf('<a href="http://third.ucllnl.org/cgi-bin/firstcutout?RA=%10.7f&Dec=%10.6f&ImageSize=1&MaxInt=1&ImageType=GIF" target="_blank">FIRST</a>',RA(I)./15,Dec(I));
   URL.FIRST{I} = sprintf('http://third.ucllnl.org/cgi-bin/firstcutout?RA=%10.7f&Dec=%10.6f&ImageSize=1&MaxInt=1&ImageType=GIF',RA(I)./15,Dec(I));

   % Link to NVSS
   Link.NVSS{I} = sprintf('<a href="http://www.cv.nrao.edu/cgi-bin/postage.pl?RA=%10.7f&Dec=%10.6f&Size=0.083333+0.083333&Equinox=2000&Type=application/postscript" target="_blank">NVSS</a>',RA(I)./15,Dec(I));
   URL.NVSS{I} = sprintf('http://www.cv.nrao.edu/cgi-bin/postage.pl?RA=%10.7f&Dec=%10.6f&Size=0.083333+0.083333&Equinox=2000&Type=application/postscript',RA(I)./15,Dec(I));
 
   % Link to DSS
   LinkDSS = get_dss(RA(I)./RAD,Dec(I)./RAD,[DSS_FOV DSS_FOV],'n','2R');
   Link.DSS2R{I} = sprintf('<a href="%s">DSS2R</a>',LinkDSS{1});
   URL.DSS2R{I}  = LinkDSS{1};

   LinkDSS = get_dss(RA(I)./RAD,Dec(I)./RAD,[DSS_FOV DSS_FOV],'n','2B');
   Link.DSS2B{I} = sprintf('<a href="%s">DSS2B</a>',LinkDSS{1});
   URL.DSS2B{I}  = LinkDSS{1};

   LinkDSS = get_dss(RA(I)./RAD,Dec(I)./RAD,[DSS_FOV DSS_FOV],'n','2I');
   Link.DSS2I{I} = sprintf('<a href="%s">DSS2I</a>',LinkDSS{1});
   URL.DSS2I{I}  = LinkDSS{1};

   LinkDSS = get_dss(RA(I)./RAD,Dec(I)./RAD,[DSS_FOV DSS_FOV],'n','1R');
   Link.DSS1R{I} = sprintf('<a href="%s">DSS1R</a>',LinkDSS{1});
   URL.DSS1R{I}  = LinkDSS{1};

   LinkDSS = get_dss(RA(I)./RAD,Dec(I)./RAD,[DSS_FOV DSS_FOV],'n','1B');
   Link.DSS1B{I} = sprintf('<a href="%s">DSS1B</a>',LinkDSS{1});
   URL.DSS1B{I}  = LinkDSS{1};

   Nrow = size(DataCell,2);
   for J=1:1:Nrow,
      CellTable{I,J} = DataCell{I,J};
   end
   CellTable{I,Nrow+1}  = Link.SDSSc{I};
   CellTable{I,Nrow+2}  = Link.SDSSo{I};
   CellTable{I,Nrow+3}  = Link.NED{I};
   CellTable{I,Nrow+4}  = Link.FIRST{I};
   CellTable{I,Nrow+5}  = Link.NVSS{I};
   CellTable{I,Nrow+6}  = Link.DSS2R{I};
   CellTable{I,Nrow+7}  = Link.DSS2B{I};
   CellTable{I,Nrow+8}  = Link.DSS2I{I};
   CellTable{I,Nrow+9}  = Link.DSS1R{I};
   CellTable{I,Nrow+10} = Link.DSS1B{I};
end


if (isempty(OutputFile)==1),
   % do not create an html file
else
   % create an html file
   html_table(OutputFile,'TableCell',CellTable,varargin{:});
end


