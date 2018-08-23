function plot_smap(Cat,varargin)
% Given a star catalog plot star map
% Package: celestial.map
% Description: Given a star catalog plot star map with optional
%              magnitude/color/type/proper motion ('/cy).
% Input  : - Catalog to plot
%          * arbitrary number of keyword and properties, where keywods are:
%            'FOV'      - Field of view [fovRA, fovDec], in [arcmin],
%                         default is field of view.
%            'Center'   - Field center [RA, Dec], default is mid coordinates.
%            'ColRA'    - Catalog column number containing RA,
%                         if three columns are given then assumed to
%                         be [H M S], otherwise [Radians].
%                         Default is 1.
%            'ColDec'   - Catalog column number containing Dec,
%                         if four columns are given then assumed to
%                         be [Sign D M S], otherwise [Radians].
%                         Default is 2.
%            'ColMag'   - Catalog column number containing magnitude.
%                         Default is 3.
%                         NaN magnitude size is set to 0.01
%            'ColMagMax'- Catalog column number containing max magnitude
%                         (e.g., variable stars), default is NaN.
%            'ColColor' - Catalog column number containing color index.
%                         If one column is given then the color index
%                         is converted to color via the 'ColorMap'
%                         property. If three columns are given then
%                         used as [r g b] color.
%                         Default is NaN (dont plot color).
%            'ColPMRA'  - Catalog column number containing PM in RA [mas/yr]
%                         Default is NaN (dont plot PM).
%            'ColPMDec' - Catalog column number containing PM in Dec [mas/yr]
%                         Default is NaN (dont plot PM).
%            'MinPM'    - Minimum proper motion [mas/yr] value to plot (as an arrow),
%                         default is NaN (don't plot PM).
%            'ColTypeE' - Catalog column number containing object edge color,
%                         default is NaN (e.g., EdgeColor is black).
%                         where object tye are:
%                         0 - EdgeColor is the same as object color.
%                         1 - EdgeColor is black.
%                         2 - EdgeColor is white.
%                         If a string with one of the above numbers is
%                         given, then set the symbol to the one given
%                         by this number.
%            'ColTypeS' - Catalog column number containing object symbol type,
%                         default is NaN (e.g., filled circles).
%                         where symbols type are:
%                         0  - filled circles
%                         10 - empty circles.
%                         11 - filled ellipse
%                         1  - empty ellipse
%                         2  - filled box
%                         12 - empty box
%                         3  - filled triangle
%                         13 - empty triangle
%                         4  - +
%                         14 - + inside colored circle
%                         24 - + inside white circle
%                         5  - X
%                         15 - X inside colored circle
%                         25 - X inside white circle
%                         94 - dotted line
%                         If a string with one of the above numbers is
%                         given, then set the symbol to the one given
%                         by this number.
%            'ColEll_A' - Catalog column number containing semi major axis of
%                         ellipse [radians], default is NaN.
%            'ColEll_B' - Catalog column number containing semi minor axis of
%                         ellipse [radians], default is NaN.
%            'ColEll_E' - Catalog column number containing ellipticity of
%                         ellipse [radians], default is NaN.
%            'ColEll_PA'- Catalog column number containing PA of
%                         ellipse [radians], default is NaN.
%            'Labels'   - cell array containing labels (one per object),
%                         default is NaN. 
%            'MagRange' - Magnitude range to map into the 'MagSize' property.
%                         Default is [min(Mag) max(Mag)].
%            'MagSize'  - Star sizes to map from star magnitude.
%                         Default is [3 30].
%            'MagFun'   - Star magnitude-size mapping table, if given 
%                         overrid MagSize.
%                         Two column matrix [Mag, Size] and use linear
%                         interpolation in between.
%                         Default is NaN.
%            'ColorMap' - Color map [plot_color, min_color, max_color],
%                         default is : [0 1 1 -Inf 0.00
%                                       0 0 1 0.00 0.75
%                                       0 1 0 0.75 1.50
%                                       1 1 0 1.50 2.25
%                                       1 0 0 2.25 3.00
%                                       1 0 1 3.00 +Inf
%                                       0 0 0 NaN  NaN
%                                       1 1 1 Inf  Inf]
%            'PlotTitle'- Plot figure title {'yes' | 'no'}, default is 'yes'.
%            'Units'    - labels units : {'rad' | 'deg' | 'arcmin' | 'arcsec' | 'Natural'},
%                         default is 'arcmin'.
%            'Project ' - map projection : {'natural'  | 'stereographic' | 'gnomonic'    |
%                                           'mercator' | 'cassini'       | 'conic'       |
%                                           'albers'   | 'bonne'         | 'cylindrical' |
%                                           'parabolic'| 'polar'         | 'sinusoidal'  |
%                                           'aitoff'   | 'hammer'        | 'mollweide'   |
%                                           'stereographic_polar' | 'azimuthal'}
% Output : null
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Mar 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%---------------------------------------------------------------------------
import celestial.coo.*
import celestial.proj.*


RAD             = 180./pi;
EllipseWidth    = 1;
MagNaNSize      = 0.01;
MagInterpMethod = 'linear';
MaxMagOuterCircColor  = 'k';
MaxMagAnnulusColor    = 'w';

Nv       = length(varargin);
if (Nv./2~=floor(Nv./2))
   error('Illigal number of unput arguments');
end

%--- set default values ---
FOV       = NaN;
Center    = NaN;
ColRA     = 1;
ColDec    = 2;
ColMag    = 3;
ColMagMax = NaN;
ColColor  = NaN;
ColPMRA   = NaN;
ColPMDec  = NaN;
MinPM     = NaN;
ColTypeE  = NaN;
ColTypeS  = NaN;
ColEll_A  = NaN;
ColEll_B  = NaN;
ColEll_E  = NaN;
ColEll_PA = NaN;
Labels    = NaN;
MagRange  = NaN;
MagSize   = [3 30];
MagFun    = NaN;
ColorMap  = [0 1 1 -Inf 0.00
             0 0 1 0.00 0.75
             0 1 0 0.75 1.50
             1 1 0 1.50 2.25
             1 0 0 2.25 3.00
             1 0 1 3.00 +Inf
             0 0 0 NaN  NaN
             1 1 1 Inf  Inf];
PlotTitle = 'yes';
Units     = 'arcmin';
Project   = 'natural';

for I=1:2:Nv-1
   switch varargin{I}
    case 'FOV'
       FOV         = varargin{I+1};
    case 'Center'
       Center      = varargin{I+1};
    case 'ColRA'
       ColRA       = varargin{I+1};
    case 'ColDec'
       ColDec      = varargin{I+1};
    case 'ColMag'
       ColMag      = varargin{I+1};
    case 'ColMagMax'
       ColMagMax   = varargin{I+1};
    case 'ColColor'
       ColColor    = varargin{I+1};
    case 'ColPMRA'
       ColPMRA     = varargin{I+1};
    case 'ColPMDec'
       ColPMDec    = varargin{I+1};
    case 'MinPM'
       MinPM       = varargin{I+1};
    case 'ColTypeE'
       ColTypeE    = varargin{I+1};
    case 'ColTypeS'
       ColTypeS    = varargin{I+1};
    case 'ColEll_A'
       ColEll_A    = varargin{I+1};
    case 'ColEll_B'
       ColEll_B    = varargin{I+1};
    case 'ColEll_E'
       ColEll_E    = varargin{I+1};
    case 'ColEll_PA'
       ColEll_PA   = varargin{I+1};
    case 'Labels'
       Labels      = varargin{I+1};
    case 'MagRange'
       MagRange    = varargin{I+1};
    case 'MagSize'
       MagSize     = varargin{I+1};
    case 'MagFun'
       MagFun      = varargin{I+1};
    case 'ColorMap'
       ColorMap    = varargin{I+1};
    case 'PlotTitle'
       PlotTitle   = varargin{I+1};
    case 'Units'
       Units       = varargin{I+1};
    case 'Project'
       Project     = varargin{I+1};
    otherwise
       error('Unknown keyword option');
   end
end

if (length(ColRA)==3)
   RA  = convertdms(Cat(:,ColRA),'H','r');
else
   RA  = Cat(:,ColRA);
end
if (length(ColDec)==4)
   Dec = convertdms(Cat(:,ColDec),'D','R');
else
   Dec  = Cat(:,ColDec);
end

MaxRA  = max(Cat(:,ColRA));
MinRA  = min(Cat(:,ColRA));
MaxDec = max(Cat(:,ColDec));
MinDec = min(Cat(:,ColDec));
if (isnan(Center)==1)
   %--- set default Center ---
   Center = [0.5.*(MaxRA+MinRA), 0.5.*(MaxDec+MinDec)];
end
if (isnan(FOV)==1)
   %--- set default Center ---
   FOV    = [MaxRA-MinRA, MaxDec-MinDec];
   %--- convert to arcmin ---
   FOV    = FOV.*RAD.*60;
else
   Dist = sphere_dist(Center(1), Center(2), RA, Dec);
   I    = find(Dist< sqrt(sum((FOV./60./RAD).^2)).*1.1);
   Cat  = Cat(I,:);
   RA   = RA(I);
   Dec  = Dec(I);
end
if (length(FOV)==1)
   FOV = [FOV, FOV];
end


if (isnan(MagRange)==1)
   %--- set default MagRange ---
   MagRange = [min(Cat(:,ColMag)), max(Cat(:,ColMag))];
end


N  = size(Cat,1);


%FieldRangeAM(1) =   - 0.5.*FOV(1)./cos(Center(2));
%FieldRangeAM(2) =   + 0.5.*FOV(1)./cos(Center(2));
%FieldRangeAM(3) =   - 0.5.*FOV(2);
%FieldRangeAM(4) =   + 0.5.*FOV(2);

%-------------------------------------------------------------------------------------------------
%--- Sort catalog by star magnitude (in order for faint star to be overplotted on bright star) ---
%-------------------------------------------------------------------------------------------------
if (isnan(ColMagMax)==1)
   Cat        = sortrows(Cat,ColMag);
else
   Cat        = sortrows(Cat,ColMagMax);
end

%----------------------
%--- set star sizes ---
%----------------------
if (numel(MagFun)==1)
   %--- Use MagSize and MagRange to set size-magnitude relation ---
   Slope      = (MagSize(1)-MagSize(2))./(MagRange(2)-MagRange(1));
   CatSize    = MagSize(1) + (Cat(:,ColMag) - MagRange(2)).*Slope;
   I          = find(CatSize<MagSize(1));
   CatSize(I) = MagSize(1);
   I          = find(CatSize>MagSize(2));
   CatSize(I) = MagSize(2);
   I          = find(isnan(CatSize)==1);
   CatSize(I) = MagNaNSize;

   if (isnan(ColMagMax)==0)
      CatSizeMax    = MagSize(1) + (Cat(:,ColMagMax) - MagRange(2)).*Slope;
      I             = find(CatSize<MagSize(1));
      CatSizeMax(I) = MagSize(1);
      I             = find(CatSizeMax>MagSize(2));
      CatSizeMax(I) = MagSize(2);
      I             = find(isnan(CatSizeMax)==1);
      CatSizeMax(I) = MagNaNSize;
   end
else
   %--- Use MagFun to set size-magnitude relation ---
   CatSize    = interp1(MagFun(:,1),MagFun(:,2),Cat(:,ColMag),MagInterpMethod);
   I          = find(isnan(CatSize)==1);
   CatSize(I) = MagNaNSize;

   if (isnan(ColMagMax)==0)
      CatSizeMax = interp1(MagFun(:,1),MagFun(:,2),Cat(:,ColMagMax),MagInterpMethod);
      I             = find(isnan(CatSizeMax)==1);
      CatSizeMax(I) = MagNaNSize;
   end
end



%----------------------
%--- set star color ---
%----------------------
if (isnan(ColColor)==1)
   %--- set color to black ---
   CatColor = zeros(N,3);
else
   %--- set star color according to its color index ---
   Ncolor = size(ColorMap,1);

   for Ic=1:1:Ncolor
      if (isnan(ColorMap(Ic,4))==1 && isnan(ColorMap(Ic,5))==1)
         I             = find(isnan(Cat(:,ColColor))==1);
         if (isempty(I)==0)
            CatColor(I,:) = ones(length(I),1)*ColorMap(Ic,[1 2 3]);
         end
      elseif (isinf(ColorMap(Ic,4))==1 && isinf(ColorMap(Ic,5))==1)
         I             = find(isinf(Cat(:,ColColor))==1);
         if (isempty(I)==0)
            CatColor(I,:) = ones(length(I),1)*ColorMap(Ic,[1 2 3]);
         end
      else
         I             = find(Cat(:,ColColor)>ColorMap(Ic,4) & Cat(:,ColColor)<=ColorMap(Ic,5));
         if (isempty(I)==0)
            CatColor(I,:) = ones(length(I),1)*ColorMap(Ic,[1 2 3]);
         end
      end
   end
end

%------------------------
%--- Ellipse settings ---
%------------------------
if (isnan(ColEll_A)==0)
   MajorAxis = Cat(:,ColEll_A);
   if (isnan(ColEll_B)==0)
      MinorAxis = Cat(:,ColEll_B);
   else
      MinorAxis = sqrt(MajorAxis.^2.*(1 - Cat(:,ColEll_E).^2));
   end
end

%-------------------
%--- PM settings ---
%-------------------
if (isnan(ColPMRA)==0 && isnan(ColPMDec)==0 && isnan(MinPM)==0)
   PlotPM   = 'y';
   CatPM    = sqrt(Cat(:,ColPMRA).^2 + Cat(:,ColPMDec).^2);
   %--- convert PM to arcmin per 100 years
   CatPMRA  = Cat(:,ColPMRA) ./(10.*60);
   CatPMDec = Cat(:,ColPMDec)./(10.*60);
else
   PlotPM   = 'n';
end

%-----------------------
%--- Marker settings ---
%-----------------------
if (isstr(ColTypeE)==1)
   %--- Set all symbols ---
   MarkerE = str2num(ColTypeE).*ones(N,1);
else
   if (isnan(ColTypeE)==1)
      MarkerE = ones(N,1);
   else
      MarkerE = Cat(:,ColTypeE);
   end
end

if (isstr(ColTypeS)==1)
   %--- Set all symbols ---
   MarkerS = str2num(ColTypeS).*ones(N,1);
else
   if (isnan(ColTypeS)==1)
      MarkerS = zeros(N,1);
   else
      MarkerS = Cat(:,ColTypeS);
   end
end

%--------------------------
%--- Plot finding chart ---
%--------------------------
figure(gcf);
%axis(FieldRangeAM);
box on;
hold on;
SexRA    = convertdms(Center(1),'r','H');
SexDec   = convertdms(Center(2),'R','D');
if (SexDec(1)==1)
   SexSign = '+';
else
   SexSign = '-';
end
switch PlotTitle
 case 'yes'
    TitleStr = sprintf('Field Center : %02d:%02d:%05.2f %c%02d:%02d:%04.1f (J2000.0)',SexRA,SexSign,SexDec(2:4)); 
    title(TitleStr);
 case 'no'
    % do nothing
 otherwise
    error('Unknown PlotTitle option');
end


%--- Select projection ---
switch Project
 case 'natural'
    % use original data - no projection
    CatRA  = Cat(:,ColRA);
    CatDec = Cat(:,ColDec);
 case 'stereographic'
    [CatRA, CatDec] = pr_stereographic(Cat(:,ColRA),Cat(:,ColDec),1,Center);
 case 'stereographic_polar'
    [CatRA, CatDec] = pr_stereographic_polar(Cat(:,ColRA),Cat(:,ColDec));
 case 'azimuthal'
    [CatRA, CatDec] = pr_azimuthal_equidist(Cat(:,ColRA),Cat(:,ColDec));
 case 'gnomonic'
    [CatRA, CatDec] = pr_gnomonic(Cat(:,ColRA),Cat(:,ColDec),1,Center);
 case 'mercator'
    [CatRA, CatDec] = pr_mercator(Cat(:,ColRA),Cat(:,ColDec),1,Center);
 case 'cassini'
    [CatRA, CatDec] = pr_cassini(Cat(:,ColRA),Cat(:,ColDec),1,Center(1));
 case 'conic'
    [CatRA, CatDec] = pr_conic(Cat(:,ColRA),Cat(:,ColDec),1,2);
 case 'albers'  
    [CatRA, CatDec] = pr_albers(Cat(:,ColRA),Cat(:,ColDec),1,Center);
 case 'bonne'
    [CatRA, CatDec] = pr_bonne(Cat(:,ColRA),Cat(:,ColDec),1,Center);
 case 'cylindrical'
    [CatRA, CatDec] = pr_cylindrical(Cat(:,ColRA),Cat(:,ColDec),1,Center);
 case 'parabolic'
    [CatRA, CatDec] = pr_parabolic(Cat(:,ColRA),Cat(:,ColDec),1);
 case 'polar'
    [CatRA, CatDec] = pr_polar(Cat(:,ColRA),Cat(:,ColDec),1);
 case 'sinusoidal'
    [CatRA, CatDec] = pr_sinusoidal(Cat(:,ColRA),Cat(:,ColDec),1);
 case 'aitoff'
    [CatRA, CatDec] = pr_aitoff(Cat(:,ColRA),Cat(:,ColDec),1);
 case 'hammer'
    [CatRA, CatDec] = pr_hammer(Cat(:,ColRA),Cat(:,ColDec),1);
 case 'mollweide'
    [CatRA, CatDec] = pr_mollweide(Cat(:,ColRA),Cat(:,ColDec),1);
 otherwise
    error('Unknown Project option');
end



%--- Select work units ---
set(gca,'XDir','reverse');
switch Project
 case 'natural'
    %--- Natural projection ---
     switch Units
      case 'arcmin'
         %---  Work in arcmin units ---
         CatRA    = (CatRA  - Center(1)).*RAD.*60;
         CatDec   = (CatDec - Center(2)).*RAD.*60;
         Hxl      = xlabel('R.A. [arcmin/cos(\delta)]');
         Hyl      = ylabel('Dec. [arcmin]');
      case 'rad'
         %---  Work in radians units ---
         CatRA    = (CatRA  - Center(1));
         CatDec   = (CatDec - Center(2));
         Hxl      = xlabel('R.A. [radians/cos(\delta)]');
         Hyl      = ylabel('Dec. [radians]');
      case 'arcsec'
         %---  Work in arcsec units ---
         CatRA    = (CatRA  - Center(1)).*RAD.*3600;
         CatDec   = (CatDec - Center(2)).*RAD.*3600;
         Hxl      = xlabel('R.A. [arcsec/cos(\delta)]');
         Hyl      = ylabel('Dec. [arcsec]');
      case 'deg'
         %---  Work in degrees units ---
         CatRA    = (CatRA  - Center(1)).*RAD;
         CatDec   = (CatDec - Center(2)).*RAD;
         Hxl      = xlabel('R.A. [deg/cos(\delta)]');
         Hyl      = ylabel('Dec. [deg]');
      case 'Natural'
         %---  Work in input natural units ---
         CatRA    = CatRA;
         CatDec   = CatDec;
         Hxl      = xlabel('');
         Hyl      = ylabel('');
      otherwise
         error('Unknown Units option');
     end
 otherwise
    %--- non natural projection ---
    %error('Not yet available');
end



%------------------
%--- Plot stars ---
%------------------
for I=1:1:N
   %--- plot marker ---
   switch MarkerS(I)
    case NaN
       if (isnan(ColMagMax)==0)
          %--- plot bigger black circle around star (e.g., variable star) ---
          Hv = plot(CatRA(I),CatDec(I),'o');
          set(Hv,'MarkerSize',CatSizeMax(I),'MarkerEdgeColor',MaxMagOuterCircColor,'MarkerFaceColor',MaxMagAnnulusColor);
       end

       H = plot(CatRA(I),CatDec(I),'o');
       set(H,'MarkerSize',CatSize(I),'MarkerFaceColor',CatColor(I,:));

    case 0
       if (isnan(ColMagMax)==0)
          %--- plot bigger black circle around star (e.g., variable star) ---
          Hv = plot(CatRA(I),CatDec(I),'o');
          set(Hv,'MarkerSize',CatSizeMax(I),'MarkerEdgeColor',MaxMagOuterCircColor,'MarkerFaceColor',MaxMagAnnulusColor);
       end

       H = plot(CatRA(I),CatDec(I),'o');
       set(H,'MarkerSize',CatSize(I),'MarkerFaceColor',CatColor(I,:));

    case 10
       if (isnan(ColMagMax)==0)
          %--- plot bigger black circle around star (e.g., variable star) ---
          Hv = plot(CatRA(I),CatDec(I),'o');
          set(Hv,'MarkerSize',CatSizeMax(I),'MarkerEdgeColor',MaxMagOuterCircColor,'MarkerFaceColor',MaxMagAnnulusColor);
       end

       H = plot(CatRA(I),CatDec(I),'o');
       set(H,'MarkerSize',CatSize(I));

    case 1  
       H = plot_ellipse([CatRA(I),CatDec(I)],[MajorAxis(I),MinorAxis(I)],[],-Cat(I,ColEll_PA)+0.5.*pi,CatColor(I,:),EllipseWidth,1./cos(Cat(I,ColDec)),[1 1 1]);
    case 11
       H = plot_ellipse([CatRA(I),CatDec(I)],[MajorAxis(I),MinorAxis(I)],[],-Cat(I,ColEll_PA)+0.5.*pi,CatColor(I,:),EllipseWidth,1./cos(Cat(I,ColDec)),CatColor(I,:));

    case 2
       H = plot(CatRA(I),CatDec(I),'s');
       set(H,'MarkerSize',CatSize(I),'MarkerFaceColor',CatColor(I,:));
    case 12
       H = plot(CatRA(I),CatDec(I),'s');
       set(H,'MarkerSize',CatSize(I));
    case 3
       H = plot(CatRA(I),CatDec(I),'^');
       set(H,'MarkerSize',CatSize(I),'MarkerFaceColor',CatColor(I,:));
    case 13
       H = plot(CatRA(I),CatDec(I),'^');
       set(H,'MarkerSize',CatSize(I));
    case 4
       H = plot(CatRA(I),CatDec(I),'k+');
       set(H,'MarkerSize',CatSize(I));
    case 14
       H = plot(CatRA(I),CatDec(I),'o');
       set(H,'MarkerSize',CatSize(I),'MarkerFaceColor',CatColor(I,:));
       H = plot(CatRA(I),CatDec(I),'k+');
       set(H,'MarkerSize',CatSize(I));
    case 24
       H = plot(CatRA(I),CatDec(I),'o');
       set(H,'MarkerSize',CatSize(I));
       H = plot(CatRA(I),CatDec(I),'k+');
       set(H,'MarkerSize',CatSize(I));
    case 5
       H = plot(CatRA(I),CatDec(I),'kx');
       set(H,'MarkerSize',CatSize(I));
    case 15
       H = plot(CatRA(I),CatDec(I),'o');
       set(H,'MarkerSize',CatSize(I),'MarkerFaceColor',CatColor(I,:));
       H = plot(CatRA(I),CatDec(I),'kx');
       set(H,'MarkerSize',CatSize(I));
    case 25
       H = plot(CatRA(I),CatDec(I),'o');
       set(H,'MarkerSize',CatSize(I));
       H = plot(CatRA(I),CatDec(I),'kx');
       set(H,'MarkerSize',CatSize(I));
    case 94
       H = plot(CatRA(I),CatDec(I),':');
       %set(H,'MarkerSize',CatSize(I));
    otherwise
       error('Unknown ColTypeS option');
   end
   %--- plot marker edge ---
    switch MarkerE(I)
    %case NaN
    %   set(H,'MarkerEdgeColor','k');
    case 0
       set(H,'MarkerEdgeColor',CatColor(I,:));
    case 1
       set(H,'MarkerEdgeColor',[0 0 0]);
    case 2
       set(H,'MarkerEdgeColor',[1 1 1]);
    otherwise
       error('Unknown ColTypeE option');
   end
   %--- plot proper motion arrow ---
   switch PlotPM
    case 'n'
       % don't plot proper motions 
    case 'y'
       if (CatPM(I)>=MinPM)
          %--- Plot PM ---
          arrow([CatRA(I),CatDec(I)],[CatRA(I)+CatPMRA(I), CatDec(I)+CatPMDec(I)],'Length',10);
       end
    otherwise
       error('Unknown PlotPM option');
   end

end
