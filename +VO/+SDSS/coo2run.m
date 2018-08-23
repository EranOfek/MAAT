function [Run,MJD,Dist]=coo2run(RA,Dec,DR,Assign)
% Convert J2000 equatorial coordinates to SDSS run/rerun/col/field ID.
% Package: VO.SDSS
% Description: Given celestial equatorial coordinates, find all SDSS
%              Run/Rerun/Col/Field number ID that cover the coordinates.
% Input  : - J2000.0 RA in [rad] or [H M S] or sexagesimal string.
%          - J2000.0 Dec in [rad] or [Sign D M S] or sexagesimal string.
%          - SDSS data relaese, default is 'DR9'.
%            If empty (i.e., []), then use default.
%          - Assign and read {0 | 1} the list of fields into/from
%            the matlab workspace.
%            Default is false.
%            If you are running coo2run many times, use 1 for faster
%            execuation.
% Output : - A cell array of lists of [Run, Rerun, Col, Field]
%            that covers the given coordinates. Ecah cell for each
%            coordinate.
%          - A cell array of modified JD of frame for [u g r i z] frames,
%            one line per run. Each cell per each coordinate.
%          - A cell array od distances (in radians) of the requested
%            RA, Dec from each side of the frame polygon.
%            Each cell per each coordinate.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Mar 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Package: sdss
% Example: [Run,MJD,D]=VO.SDSS.coo2run([10 04 00],[1 41 0 0]);
%          min(D{1}.*RAD.*3600,[],2)  % print the distance ["] of the
%                                     % coordinate from nearest boundry.
% Reliable: 2
%--------------------------------------------------------------------------

CatField    = AstCat.CatField;
ColField    = AstCat.ColField;

BaseVarName = 'AstCat_SDSS_Fields';

RAD    = 180./pi;
Threshold  = 0.3./RAD;   % search threshold (larger than frame size)

DefDR = 'DR9';
Def.Assign = false;
if (nargin==2)
   DR     = DefDR;
   Assign = Def.Assign;
elseif (nargin==3)
   Assign = Def.Assign;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(DR)==1)
   DR = DefDR;
end
FileName    = sprintf('SDSS_%s_Fields_All_PolySort.hdf5',DR);



RA  = celestial.coo.convertdms(RA,'gH','r');
Dec = celestial.coo.convertdms(Dec,'gD','R');

if (Assign)
    try
       Fields = evalin('base',sprintf('%s;',BaseVarName));
    catch

       %--- Load SDSS Fields ---
       %AstCat_SDSS_Fields = AstCat.loadh2astcat(FileName);
       AstCat_SDSS_Fields = cats.SDSS.SDSS_DR9_Fields_All_PolySort;
%        eval(sprintf('load %s.mat',FieldsStr));
%        eval(sprintf('Fields = %s;',FieldsStr));   % copy to Fields
%        eval(sprintf('clear %s;',FieldsStr));

       assignin('base',BaseVarName,AstCat_SDSS_Fields);
    end
else
    %--- Load SDSS Fields ---
    %AstCat_SDSS_Fields = AstCat.loadh2astcat(FileName);
    AstCat_SDSS_Fields = cats.SDSS.SDSS_DR9_Fields_All_PolySort;

end

Col = AstCat_SDSS_Fields.(ColField);


Ind  = [1 2; 2 3; 3 4; 4 1];   % indecies of sides
N    = length(RA);
Run  = cell(N,1);
MJD  = cell(N,1);
Dist = cell(N,1);
L = search_cat(AstCat_SDSS_Fields.(CatField)(:,[10 14]),RA(:),Dec(:),'SearchRad',Threshold);
for I=1:1:N
   %L = cat_search(Fields,[ColRA ColDec],[RA(I) Dec(I)],Threshold,'circle','Dec','sphere');

   InL = [];
   for J=1:1:L(I).Nfound
      Inside1 = inpolygon(RA(I),Dec(I),AstCat_SDSS_Fields.(CatField)(L(I).IndCat(J),Col.RA1+[0:1:3]), ...
                                       AstCat_SDSS_Fields.(CatField)(L(I).IndCat(J),Col.Dec1+[0:1:3]));
      Inside2 = inpolygon(RA(I)-2.*pi,Dec(I),AstCat_SDSS_Fields.(CatField)(L(I).IndCat(J),Col.RA1+[0:1:3]), ...
                                             AstCat_SDSS_Fields.(CatField)(L(I).IndCat(J),Col.Dec1+[0:1:3]));
      Inside3 = inpolygon(RA(I)+2.*pi,Dec(I),AstCat_SDSS_Fields.(CatField)(L(I).IndCat(J),Col.RA1+[0:1:3]), ...
                                             AstCat_SDSS_Fields.(CatField)(L(I).IndCat(J),Col.Dec1+[0:1:3]));
      Inside  = sign(Inside1 + Inside2 + Inside3);
      if (Inside1 || Inside2 || Inside3)
         InL = [InL; L(I).IndCat(J)];
      end
   end
   Run{I} = AstCat_SDSS_Fields.(CatField)(InL,[Col.Run, Col.Rerun, Col.Camcol, Col.Field]);
   MJD{I} = AstCat_SDSS_Fields.(CatField)(InL,[Col.MJD_u, Col.MJD_g, Col.MJD_r, Col.MJD_i, Col.MJD_z]);

   %--- calculate distance of point from boundries ---
   % assuming plane(!) geometry (approximation)
   for Is=1:1:length(Ind)
      ColRA1  = Col.RA1  - 1 + Ind(Is,1);
      ColRA2  = Col.RA1  - 1 + Ind(Is,2);
      ColDec1 = Col.Dec1 - 1 + Ind(Is,1);
      ColDec2 = Col.Dec1 - 1 + Ind(Is,2);

      D0 = celestial.coo.sphere_dist(AstCat_SDSS_Fields.(CatField)(InL,ColRA1), ...
                       AstCat_SDSS_Fields.(CatField)(InL,ColDec1), ...
                       AstCat_SDSS_Fields.(CatField)(InL,ColRA2), ...
                       AstCat_SDSS_Fields.(CatField)(InL,ColDec2));
      D1 = celestial.coo.sphere_dist(AstCat_SDSS_Fields.(CatField)(InL,ColRA1), ...
                       AstCat_SDSS_Fields.(CatField)(InL,ColDec1), ...
                       RA(I), Dec(I));
      D2 = celestial.coo.sphere_dist(AstCat_SDSS_Fields.(CatField)(InL,ColRA2), ...
                       AstCat_SDSS_Fields.(CatField)(InL,ColDec2), ...
                       RA(I), Dec(I));

      S = 0.5.*(D0+D1+D2);    % semi-perimeter
      H = 2.*sqrt(S.*(S-D0).*(S-D1).*(S-D2))./D0;  % Height using Heron formula
      Dist{I}(:,Is) = H;
   end
end
clear Fields;

