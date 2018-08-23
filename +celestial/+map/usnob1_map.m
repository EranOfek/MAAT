function FieldRangeDeg=usnob1_map(RA,Dec,FOV,Color,Gal,Stars,MagRange,SymbolSize,GalFaceColor)
% Plot a finding chart using a local copy of the USNO-B2.0
% Package: celestial.map
% Description: Plot a finding chart using a local copy of the USNO-B2.0
%              catalog. The user can select between b/w stars or color
%              stars (with their O-E color index color coded).
%              If O-E is not available then stars are plotted in black.
%              The user can overplot known galaxies from the PGC catalog.
%              In the color option the edge of probable stellar objects
%              is marked by black circle. The user can overplot known 
%              bright stars (VT<11) for which spikes and saturation
%              artifact can be seen.
% Input  : - RA (J2000.0), [HH MM SS] or [Radians]
%          - Dec (J2000.0), [Sign, DD MM SS] or [Radians]
%          - Field of view [arcsec], default is 300 arcsec.
%          - Color option: 'bw'-black & white; 'sp'-spectral type (default)
%          - Add PGC galaxies, default is 'y'.
%          - Add known bright stars with an estimat of their artifact region.
%          - Magnitude range, default is [-2 22].
%          - Symbol plot size range [3, 30].
%          - Galaxy face color, ('w' - will delete stars found behind the galaxy),
%            default is 'none'.
% Output : - Field range in deg [RAmin, RAmax, Decmin, Decmax].
% Plot   : Finding chart.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                     March 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: apmread.m, get_usnob1sel.m, pgc.dat, tyc2_11.dat
%    Need: ubcone (binary file!) + USNO-B2.0 catalog on disk
% Reliable: NOT WORKING ANYMORE
%------------------------------------------------------------------------------
RAD   = 180./pi;

import celestial.coo.*

%--- additional options ---
PlotLegend      = 'y';
GalaxyColor     = [0 0 0];
GalaxyWidth     = 1;
SaturationColor = [1 0 0];
SaturationWidth = 2;
LegendMag       = [22; 20; 18; 16; 14];

%--- spectral type definition ---
SymColorR       =   [0.00; 0.75; 1.50; 2.25; 3.00; NaN];
SymColor        = ['c';  'b' ; 'g' ; 'y' ; 'r' ; 'm'; 'k'];
SymColorText    = [-Inf; SymColorR(1:end-1); Inf; SymColorR(end)];

%--- Default values ---
DefFOV          = 300;
DefColor        = 'sp';
DefGal          = 'y';
DefStars        = 'y';
DefMagRange     = [-2 22];
DefSymbolSize   = [3, 30];
DefGalFaceColor = 'none';
if (nargin==2),
   FOV         = DefFOV;
   Color       = DefColor;
   Gal         = DefGal;
   Stars       = DefStars;
   MagRange    = DefMagRange;
   SymbolSize  = DefSymbolSize;
   GalFaceColor= DefGalFaceColor;
elseif (nargin==3),
   Color       = DefColor;
   Gal         = DefGal;
   Stars       = DefStars;
   MagRange    = DefMagRange;
   SymbolSize  = DefSymbolSize;
   GalFaceColor= DefGalFaceColor;
elseif (nargin==4),
   Gal         = DefGal;
   Stars       = DefStars;
   MagRange    = DefMagRange;
   SymbolSize  = DefSymbolSize;
   GalFaceColor= DefGalFaceColor;
elseif (nargin==5),
   Stars       = DefStars;
   MagRange    = DefMagRange;
   SymbolSize  = DefSymbolSize;
   GalFaceColor= DefGalFaceColor;
elseif (nargin==6),
   MagRange    = DefMagRange;
   SymbolSize  = DefSymbolSize;
   GalFaceColor= DefGalFaceColor;
elseif (nargin==7),
   SymbolSize  = DefSymbolSize;
   GalFaceColor= DefGalFaceColor;
elseif (nargin==8),
   GalFaceColor= DefGalFaceColor;
elseif (nargin==9),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (length(RA)==3),
   RA  = convertdms(RA,'H','r');
end
if (length(Dec)==4),
   Dec = convertdms(Dec,'H','R');
end


Col_RA        = 1;
Col_Dec       = 2;
Col_B1        = 15;
Col_R1        = 22;
Col_B2        = 29;
Col_R2        = 36;
Col_I2        = 43;
Col_sgB1      = 18;
Col_sgR1      = 25;
Col_sgB2      = 32;
Col_sgR2      = 39;
Col_sgI2      = 46;

Cat = get_usnob1sel(RA,Dec,FOV,'b');
%--- sort by mag ---
Cat = sortrows(Cat,Col_R2);

I   = find((Cat(:,Col_R2)>=MagRange(1) & Cat(:,Col_R2)<=MagRange(2)) | ...
           (Cat(:,Col_R1)>=MagRange(1) & Cat(:,Col_R1)<=MagRange(2)));
Cat = Cat(I,:);
N   = size(Cat,1);

%FieldRangeDeg(1) = RA.*RAD  - 0.5.*FOV./(3600.*cos(Dec));
%FieldRangeDeg(2) = RA.*RAD  + 0.5.*FOV./(3600.*cos(Dec));
%FieldRangeDeg(3) = Dec.*RAD - 0.5.*FOV./3600;
%FieldRangeDeg(4) = Dec.*RAD + 0.5.*FOV./3600;
FieldRangeDeg(1) =   - 0.5.*FOV./(3600.*cos(Dec));
FieldRangeDeg(2) =   + 0.5.*FOV./(3600.*cos(Dec));
FieldRangeDeg(3) =   - 0.5.*FOV./3600;
FieldRangeDeg(4) =   + 0.5.*FOV./3600;

%--- prepare catalog ---
SG       = nanmean(Cat(:,[Col_sgB1 Col_sgR1 Col_sgB2 Col_sgR2 Col_sgI2]).')';
CatType  = zeros(N,1); % galaxies
I        = find(SG>=3);
CatType(I) = 1;    % stars

BR1      = Cat(:,Col_B1) - Cat(:,Col_R1);
BR2      = Cat(:,Col_B2) - Cat(:,Col_R2);
CatBR    = BR2;
I        = find(isnan(BR2)==1);
CatBR(I) = BR1(I);

switch Color
 case 'bw'
    CatColor = zeros(N,3);
 case 'sp'
    I1       = find(CatBR<SymColorR(1));
    I2       = find(CatBR>=SymColorR(1) & CatBR<SymColorR(2));
    I3       = find(CatBR>=SymColorR(2) & CatBR<SymColorR(3));
    I4       = find(CatBR>=SymColorR(3) & CatBR<SymColorR(4));
    I5       = find(CatBR>=SymColorR(4) & CatBR<SymColorR(5));
    I6       = find(CatBR>=SymColorR(5));
    In       = find(isnan(CatBR)==1);
    CatColor(I1,:) = SymColor(1,:);
    CatColor(I2,:) = SymColor(2,:);
    CatColor(I3,:) = SymColor(3,:);
    CatColor(I4,:) = SymColor(4,:);
    CatColor(I5,:) = SymColor(5,:);
    CatColor(I6,:) = SymColor(6,:);
    CatColor(In,:) = SymColor(7,:);
 otherwise
    error('Unknown Color option');
end


CatMag   = Cat(:,Col_R2);
I        = find(isnan(CatMag)==1);
CatMag(I)= Cat(I,Col_R1);
% stars fainter the MagRange(2) - set magnitude to MagRange(2)
I        = find(CatMag>MagRange(2));
CatMag(I)= MagRange(2);
Slope    = (SymbolSize(1)-SymbolSize(2))./(MagRange(2)-MagRange(1));
CatSize  = SymbolSize(1) + (CatMag - MagRange(2)).*Slope;

CatRA    = Cat(:,1).*RAD;
CatDec   = Cat(:,2).*RAD;

%--- Plot finding chart ---
figure;
axis(FieldRangeDeg.*60);
box on;
hold on;
SexRA    = convertdms(RA,'r','H');
SexDec   = convertdms(Dec,'R','D');
if (SexDec(1)==1),
   SexSign = '+';
else
   SexSign = '-';
end
TitleStr = sprintf('Field Center : %02d:%02d:%05.2f %c%02d:%02d:%04.1f (J2000.0)',SexRA,SexSign,SexDec(2:4)); 
title(TitleStr);
CatRA    = (CatRA - RA.*RAD).*60;
CatDec   = (CatDec - Dec.*RAD).*60;
Hxl      = xlabel('R.A. [arcmin/cos(\delta)]');
Hyl      = ylabel('Dec. [arcmin]');
set(gca,'XDir','reverse');


Hstars = zeros(N,1);
for I=1:1:N,
   H = plot(CatRA(I),CatDec(I),'o');
   Hstars(I) = H;
   switch Color
    case 'sp'
       if (CatType(I)==0),
          % galaxy
          set(Hstars(I),'MarkerSize',CatSize(I),'MarkerEdgeColor',CatColor(I),'MarkerFaceColor',CatColor(I));
       else
          % star
          set(Hstars(I),'MarkerSize',CatSize(I),'MarkerEdgeColor','k','MarkerFaceColor',CatColor(I));
       end
    case 'bw'
       set(Hstars(I),'MarkerSize',CatSize(I),'MarkerEdgeColor','k','MarkerFaceColor','k');
    otherwise
       error('Unknown Color option');
   end
end


%-------------------------
%--- overplot galaxies ---
%-------------------------
switch Gal
 case 'y'
    load('pgc.mat'); 
    PGC = pgc;
    clear pgc;
	
    IndPGC  = cat_search(PGC,[1 2],[RA Dec],sqrt(2).*FOV./(3600.*RAD));
    GalCoo    = [PGC(IndPGC,1) - RA, PGC(IndPGC,2) - Dec].*RAD.*60;
    MajorAxis = 0.5.*0.1.*10.^PGC(IndPGC,4);
    MinorAxis = MajorAxis./(10.^PGC(IndPGC,6));
    PA        = PGC(IndPGC,8)./RAD;

    Ng        = size(GalCoo,1);
    for Ig=1:1:Ng,
       H = plot_ellipse(GalCoo(Ig,1:2),[MajorAxis(Ig),MinorAxis(Ig)],[],-PA(Ig)+0.5.*pi,GalaxyColor,GalaxyWidth,cos(Dec),GalFaceColor);
    end
 case 'n'
    % do nothing
 otherwise
    error('Unknown Gal option');
end



%-----------------------------
%--- overplot bright stars ---
%-----------------------------
switch Stars
 case 'y'
    load('tyc2_11.mat'); 
    Tyc2 = tyc2_11;
    clear tyc2_11;

    IndS    = cat_search(Tyc2,[1 2],[RA Dec], 2.*sqrt(2).*FOV./(3600.*RAD));
    Ns      = length(IndS);

    StarCat = [[Tyc2(IndS,1)-RA, Tyc2(IndS,2)-Dec].*RAD.*60, Tyc2(IndS,7)];

    for Is=1:1:Ns,
       [Circ,Spike] = saturation_region(StarCat(Is,3));
       Circ         = Circ./60;
       Spike        = Spike./60;
       Hc = plot_ellipse(StarCat(Is,1:2),[Circ,Circ],[],0,SaturationColor,SaturationWidth,cos(Dec));
       SpikesX = [[StarCat(Is,1)+ [-Spike;+Spike]./cos(Dec)], StarCat(Is,2)*ones(2,1)];
       Hsx = plot(SpikesX(:,1),SpikesX(:,2));
       set(Hsx,'Color',[0 1 0],'LineWidth',2);
       SpikesY = [StarCat(Is,1).*ones(2,1), StarCat(Is,2)+[+Spike;-Spike]];
       Hsy = plot(SpikesY(:,1),SpikesY(:,2));
       set(Hsy,'Color',[0 1 0],'LineWidth',2);

    end

 case 'n'
    % do nothing
 otherwise
    error('Unknown Stars option');
end

%-------------------
%--- Plot legend ---
%-------------------
switch PlotLegend
 case 'y'
    % given the default aspect ratio plot box
    Hmap = gca;
    set(gca,'Position',[0.11 0.11 0.58 0.775]);

    Hlegend = axes('position',[0.70 0.11 0.25 0.775]);
    axis([0 1 0 1]);
    hold on;    

    %--- magnitude symbol size ---
    LegendMagSize  = SymbolSize(1) + (LegendMag - MagRange(2)).*Slope;
    Nl             = length(LegendMag);
    OffsetMagText = 0.3;
    XposL = 0.8;
    Ydec  = 0.07;
    YposL = 1.0;
    Ht = text(1.0,YposL-Ydec.*1.5,'E-band mag');
    set(Ht,'Rotation',-90);
    for Il=1:1:Nl,
       YposL = YposL - Ydec;
       Hl = plot(XposL,YposL,'o');
       set(Hl,'MarkerSize',LegendMagSize(Il),'MarkerEdgeColor','k','MarkerFaceColor','k');
       Ht = text(XposL-OffsetMagText,YposL,sprintf('%04.1f',LegendMag(Il)));
       set(Ht,'FontSize',10);
   end

   %--- color legend ---
   Nc = length(SymColor);
   YposL = YposL - 0.08;
   Ht = text(1.0,YposL-Ydec.*2,'O-E color');
   set(Ht,'Rotation',-90);
   for Ic=1:1:Nc,
      YposL = YposL - Ydec.*0.8;
      Hc = plot(XposL,YposL,'o');
      set(Hc,'MarkerSize',LegendMagSize(3),'MarkerEdgeColor',SymColor(Ic),'MarkerFaceColor',SymColor(Ic));
      if (isnan(SymColorText(Ic))==1 | isnan(SymColorText(Ic+1))==1),
         Ht = text(XposL-OffsetMagText.*1.6,YposL,sprintf('Unknown'));
      else
         Ht = text(XposL-OffsetMagText.*1.6,YposL,sprintf('%4.2f..%4.2f',SymColorText(Ic),SymColorText(Ic+1)));
      end

   end

   %--- plot galaxy ---
   YposL = YposL - Ydec;
   plot_ellipse([XposL,YposL],Ydec.*0.8,0.97,10./RAD);
   Ht = text(XposL-OffsetMagText.*1.6,YposL,sprintf('galaxy'));


   set(Hlegend,'Visible','off')

 case 'n'
    FigPos = get(gcf,'Position');
    set(gcf,'Position',[FigPos(1), FigPos(2).*0.75, FigPos(3), FigPos(4).*1.26])
 otherwise
    error('Unknown PlotLegend option');
end

%--- return focus to map gca ---
set(gcf,'CurrentAxes',Hmap);



function [Circ,Spike]=saturation_region(VT);
%--- Given star VT magnitude return saturation region and spike length in arcsec ---

Circ  = 10.^(-0.11939.*VT +       2.3211);
Spike = 10.^(-0.20582.*VT +       3.4654);
