function plot_monthly_smap(Date,varargin)
% Plot a monthly sky map
% Package: celestial.map
% Description: Plot a monthly sky map with a naked eye stars for a given
%              time and observer Geodetic position. Optionaly mark planets
%              position, constellations and the milky way.
% Input  : - Date [D M Y H M S]  | 
%            Date [D M Y Frac] |
%            JD.
%          * arbitrary number of keyword and properties, where keywods are:
%            'GeodPos'    - Geodetic poistion
%                           [East Long, North Lat, Height(m)] 
%                           in radians, default is [35 32 0]./RAD.
%            'Cat'        - Star Catalog, default is mag6.mat 
%            'ColRA'      - RA column in star catalog, default is 1.
%            'ColDec'     - RA column in star catalog, default is 2.
%            'ColMag'     - Mag column in star catalog, default is 3.
%            'ColCol'     - Color (B-V) column in star catalog,
%                           default is 4.
%            'ColPMRA'    - PM RA*cos(dec) column in star catalog,
%                           default is 5 [mas/yr].
%            'ColPMDec'   - PM Dec column in star catalog,
%                           default is 6 [mas/yr].
%            'MW'         - Plot MilkyWay {'yes' | 'no'}, default is 'yes'.
%            'StarSizeScale'-Scale factor for stars size (scales linearly
%                           'MagFun'), default is 0.7.
%            'MagFun'     - Two column matrix relating star mag (first
%                           column) to symbol size (second column).
%            'StarType    - Star symbols {'black-edge' | 'same-edge'},
%                           default is 'black-edge'.
%            'ConLines'   - Plot Constellation lines {'yes' | 'no'},
%                           default is 'yes'.
%            'ConLinesColor'- Constellation lines color,
%                           default is [0 0 1]. 
%            'ConLinesWidth'- Constellation lines width, default is 0.5.
%            'ConLabels'  - Constellation labels {'yes' | 'no'},
%                           default is 'yes'.
%            'ConLabelsColor'- Constellations labels color, default is 'k'.
%            'ConLabelsOrient'- Constellations labels orientation:
%                           {'south' | 'radial'}, default is 'south'.
%            'MagLimit'   - Mag limit, default is 5.
%            'PM'         - Applay proper motion to star catalog
%                           {'yes' | 'no'}, default is 'yes'.
%            'Precess'    - Precess coordinates, {'yes' | 'no'},
%                           default is 'yes'.
%            'Equinox'    - Star Catalog mean equinox in JD,
%                           default is 2451545.5
%            'Epoch'      - Star Catalog epoch in JD, default is 2451545.5
%            'ColorBrightMW' - Color scheme for Milky way bright patches,
%                           default is [0.8 0.8 0.8].
%            'ColorDarkMW' - Color scheme for Milky way dark patches,
%                           default is [1 1 1].
%            'Planets'    - Show planets, the following options are available:
%                           If string:
%                           'no'/'off' - don't show planets - default.
%                           'ALL' - show Mercury..Neptune + Moon & Sun.
%                           'All' - show Mercury..Neptune.
%                           'EYE' - Show Mecury..Saturn + Moon & Sun.
%                           'Eye' - Show Mecury..Saturn.
%                           Or vector of numbers, one number for each planets
%                           according to the following codes:
%                           1-Mercury; 2-Venus; 3-Moon; 4-Mars;
%                           5-Jupiter; 6-Saturn; 7-Uranus; 8-Neptune;
%                           0-Sun.
%            'PlJD'       - Vector of JD in which to show planets position
%                           (e.g., multiple epochs), default is the single
%                           epoch given in Date parameter.
%            'PlMarker'   - Planet marker type, default is 'o'.
%            'PlLine'     - Planet position connecting line style,
%                           default is '-'.
%            'PlMarkSize' - Planet marker size, default is 'Mag'.
%                          (if 'Mag' then size according to magnitude scheme
%                           of stars.)
%                           If 'Icon', then use single size planets icons.
%            'PlColor'    - Planets color, default is [0 0 0].
%            'MoonPhases' - Plot Moon phases calendar for current month,
%                           {'y','n'}, default is 'n'.
%            'ObjectList' - Matrix of additional coordinates/objects to plot:
%                           [RA, Dec, Equinox(JD), LineType, SymbolType,
%                            Mag, SymbolSize,
%                            SymbolColor(3columns), TailPA, TailLength].
%                           where different set of objects may be entered
%                           with a lines of NaNs seprating between them.
%                           LineType   : 0-no line, 1-solid, 2-dashed,
%                                        3-dash-dot, 4-dpineotted
%                           SymbolType : 0-no symbol, 1-'o', 2-'^', 3-'s',
%                                        4-'p', 5-'h',
%                                        6-comet, 7-meteor shower
%                           Mag        : plot symbol size according to
%                                        magnitude.
%                           SymbolSize : If symbol size not NaN then used
%                                        the value in
%                                        SymbolSize column instead of Mag.
%                           SymbolColor: Three column vector of symbol color
%                           TailPA     : Comet tail P.A. in radians.
%                           TailLength : Comet tail length in radians.
%            'ColorOut'   - Color outside the map, default is [1 1 1]. 
%            'ColorIn'    - Color inside the map, default is [1 1 1].
%            'ColorScheme'- Select a predefined color scheme
%                           (override user selected colors):
%                           'C'   - Colored stars, black planets,
%                                   white background,
%                                   gray milky way, default.
%                           'BW'  - Black stars and planets, white bacground, 
%                                   gray milky way.
%                           'CB'  - Color stars, white planets,
%                                   black background, 
%                                   blue milky way.
%            'Copyright'  - {'off' | 'EO' | 'TAU'}, default is 'TAU'.
%            'Legend'     - {'off' | 'Mag'}, default is 'Mag'.
%            'Visible'    - {'on' | 'off'}, default is 'off'
%                           (for plot visibility).
%            'Ecliptic'   - plot ecliptic {'yes'|'no'}, default is 'no'.
% Output : null
% Plot   : Plot a monthly sky map.
% Tested : Matlab 6.5
%     By : Eran O. Ofek                    Nov 2004
% Web example: http://astroclub.tau.ac.il/skymaps/monthly/
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.map.plot_monthly_smap([1 1 2010 18./24]);
% Reliable: 1
%--------------------------------------------------------------------------
import celestial.coo.*
import celestial.time.*
import celestial.map.*
import celestial.SolarSys.*



RAD        = 180./pi;
ARCSEC_DEG = 3600;
JulianYear = 365.25;
DeltaT     = 0;   % in days
AltLimit   = 0;
MapRadius  = 1;

LegendStarColor = [0 1 0];

import celestial.proj.*
import celestial.stars.*
import celestial.time.*
import celestial.SolarSys.*
import celestial.coo.*
import celestial.Earth.*

% M2
MagFun   = [-30 12;...
             -5 7;...
              0 6;...
              1 5;...
              2 4;...
              3 3;...
              4 2;...
              5 1;...
              9 1];


% M3
MagFun   = [-30 12;...
             -5 9;...
              0 8;...
              1 7;...
              2 5.5;...
              3 4;...
              4 2.5;...
              5 1;...
              9 1];

% M4
MagFun   = [-30 12;...
             -5 8;...
              0 6.5;...
              1 5;...
              2 4;...
              3 3;...
              4 2;...
              5 1;...
              9 1];

% M1
MagFun   = [-30 12;...
             -5 12;...
              0 10;...
              1 9;...
              2 7;...
              3 5;...
              4 4;...
              5 1;...
              9 1];


% M5
MagFun   = [-30 12;...
             -5 9;...
              0 7.5;...
              1 5.5;...
              2 4;...
              3 3;...
              4 2;...
              5 1;...
              9 1];

% M6
MagFun   = [-30 12;...
             -5 10;...
              0 8.5;...
              1 6.5;...
              2 4;...
              3 3;...
              4 2;...
              5 1;...
              9 1];


%--- old color map ---
%ColorMapC = [0 0 1   -Inf 0.20;...
%             0 1 1   0.20 0.60;...
%             0 1 0   0.60 1.00;...
%             1 1 0   1.00 1.30;...
%             1 0 0   1.30 1.60;...
%             1 0 1   1.60 5.00;...
%             0 0 0   NaN NaN;...
%             1 1 1   Inf Inf];

% New color map corresponding to B-V color index of...
% O B A F G K M
ColorMapC = [0 0   1   -Inf  -0.30;...
             0 0.5 1   -0.30 +0.00;...
             0 1   1   +0.00 +0.37;...
             0 1   0   +0.37 +0.70;...
             1 1   0   +0.70 +1.10;...
             1 0.5 0   +1.10 +1.70;...
             1 0   0   +1.70 +Inf;...
             0 0   0    NaN   NaN];

ColorMapB = [0 0 0 -Inf -10;...
             0 0 0 -10  +20;...
             0 0 0 +20  +Inf];


%M6 - cluster
%M7 - cluster
%h&c Per - cluster
%M31 - gal
%M33 - gal
%M42 - neb
%M13 - glob
%M45 - cluster
%M8  - neb
%M11 - cluster
%M16 - nebula
%M17 - nebula
%M44 - cluster
%M22 - glob
%Omega Cen - glob

NebulaeCat = [       265.07        -32.2
       268.45      -34.783
       34.517       57.516
       10.683       41.269
       23.463        30.66
       83.821      -5.3911
       250.42        36.46
        56.85       24.117
        270.9      -24.387
       282.75      -6.2667
        274.7      -13.817
       275.11      -16.177
        130.1       19.683
        279.1      -23.903
	201.69      -47.477]./RAD;
Lnc=size(NebulaeCat,1);
NebulaeCat = [NebulaeCat, ones(Lnc,1).*3, zeros(Lnc,11)];

%NebulaeCat(1,:) = [[10.6750, 41.2667]./RAD 3 0 zeros(1,10)];   % M31
%NebulaeCat(2,:) = [[83.8221, -5.3911]./RAD 3 0 zeros(1,10)];   % M42
%NebulaeCat(3,:) = [[23.4750, 30.6500]./RAD 3 0 zeros(1,10)];   % M33
%NebulaeCat(4,:) = [[201.6912, -47.4769]./RAD 3 0 zeros(1,10)]; % omega Cen
%NebulaeCat(5,:) = [[250.4227, 36.4603]./RAD 3 0 zeros(1,10)];  % M13
%NebulaeCat(6,:) = [[270.8625, -24.3867]./RAD 3 0 zeros(1,10)]; % M8
%NebulaeCat(7,:) = [[282.7500, -6.2667]./RAD 3 0 zeros(1,10)];  % M11
%NebulaeCat(8,:) = [[130.1000, 19.6833]./RAD 3 0 zeros(1,10)];  % M44
%NebulaeCat(9,:) = [[279.1009, -23.9034]./RAD 3 0 zeros(1,10)]; % M22
%NebualeCat(10,:)= [[35.1375, 57.1417]./RAD 3 0 zeros(1,10)];   % h & chi Per
% M45 - alreday marked


%-----------------------------------
%--- Read date and convert to JD ---
%-----------------------------------
DateL = length(Date);
if (DateL==1)
   JD = Date;
elseif (DateL==4)
   JD = julday(Date);
elseif (DateL==6)
   Frac = convertdms(Date(:,4:6),'H','f');
   JD   = julday([Date(:,1:3), Frac]);
else
   error('Unkown Date Format');
end

%--------------------------------
%--- Read varargin parameters ---
%--------------------------------
Nv       = length(varargin);
if (Nv./2~=floor(Nv./2))
   error('Illigal number of unput arguments');
end

%--- set default values ---
GeodPos     = [35 32 0]./RAD;
Cat         = 'cats.bright.mag6';
ColRA       = 1;
ColDec      = 2;
ColMag      = 3;
ColCol      = 4;
ColPMRA     = 5;
ColPMDec    = 6;
MW          = 'yes';
StarType    = 'black-edge';
StarSizeScale = 0.7;
ConLines    = 'yes';
ConLinesColor = [0 0 1];
ConLabels   = 'yes';
ConLabelsColor = 'k';
ConLinesWidth  = 0.5;
ConLabelsOrient = 'south';
MagLimit    = 5;
PM          = 'yes';
Precess     = 'yes';
Equinox     = 2451545.5;
Epoch       = 2451545.5;
ColorBrightMW = [0.8 0.8 0.8];
ColorDarkMW = [1 1 1];
Planets     = 'no';
PlJD        = JD;
PlMarker    = 'o';
PlLine      = '-';
PlMarkSize  = 'Mag';
PlColor     = [0 0 0];
MoonPhases  = 'n';
ObjectList  = [];
ColorOut    = [1 1 1];
ColorIn     = [1 1 1];
ColorScheme = 'C';
Copyright   = 'TAU';
Legend      = 'Mag';
Visible     = 'on';
Ecliptic    = 'no';

for I=1:2:Nv-1
   switch lower(varargin{I})
    case 'geodpos'
       GeodPos        = varargin{I+1};
    case 'cat'  
       Cat            = varargin{I+1};
    case 'colra'
       ColRA          = varargin{I+1};
    case 'coldec'
       ColDec         = varargin{I+1};
    case 'colmag'
       ColMag         = varargin{I+1};
    case 'colcol'
       ColCol         = varargin{I+1};
    case 'colpmra'
       ColPMRA        = varargin{I+1};
    case 'colpmdec'
       ColPMDec       = varargin{I+1};
    case 'mw'
       MW             = varargin{I+1};
    case 'startype'
       StarType       = varargin{I+1};
    case 'starsizescale'
       StarSizeScale  = varargin{I+1};
    case 'magfun'
       MagFun         = varargin{I+1};
    case 'conlines'
       ConLines       = varargin{I+1};
    case 'conlinescolor'
       ConLinesColor  = varargin{I+1};
    case 'conlineswidth'
       ConLinesWidth  = varargin{I+1};
    case 'conlabels'
       ConLabels      = varargin{I+1};
    case 'conlabelscolor'
       ConLabelsColor = varargin{I+1};
    case 'conlabelsorient'
       ConLabelsOrient= varargin{I+1};
    case 'maglimit'
       MagLimit       = varargin{I+1};
    case 'pm'
       PM             = varargin{I+1};
    case 'precess'
       Precess        = varargin{I+1};
    case 'equinox'
       Equinox        = varargin{I+1};
    case 'epoch'
       Epoch          = varargin{I+1};
    case 'colorbrightmw'
       ColorBrightMW  = varargin{I+1};
    case 'colordarkmw'
       ColorDarkMW    = varargin{I+1};
    case 'planets'
       Planets        = varargin{I+1};
    case 'pljd'
       PlJD           = varargin{I+1};
    case 'plmarker'
       PlMarker       = varargin{I+1};
    case 'plline'
       PlLine         = varargin{I+1};
    case 'plmarksize'
       PlMarkSize     = varargin{I+1};
    case 'plcolor'
       PlColor        = varargin{I+1};
    case 'moonphases'
       MoonPhases      = varargin{I+1};
    case 'objectlist'
       ObjectList     = varargin{I+1};
    case 'colorOut'
       ColorOut       = varargin{I+1};
    case 'colorIn'
       ColorIn        = varargin{I+1};
    case 'colorScheme'
       ColorScheme    = varargin{I+1};
    case 'copyright'
       Copyright      = varargin{I+1};
    case 'legend'
       Legend         = varargin{I+1};
    case 'visible'
       Visible        = varargin{I+1};
    case 'ecliptic'
       Ecliptic        = varargin{I+1};

    otherwise
       error('Unknown keyword option');
   end
end

%---------------------------
%--- Scaling Stars Sizes ---
%---------------------------
MagFun(:,2) = MagFun(:,2).*StarSizeScale;

%------------------------
%--- Set Color Scheme ---
%------------------------
switch ColorScheme
 case 'C'
    ColorMap        = ColorMapC;
    ColorBrightMW   = [0.8 0.8 0.8];    % gray
    ColorDarkMW     = [1 1 1];          % white
    PlColor         = [0 0 0];          % black
    ConLinesColor   = [0 0 1];          % blue
    ConLabelsColor  = [0 0 0];          % black
    ColorIn         = [1 1 1];          % white 
    ColorOut        = [1 1 1];          % white
 case 'BW'
    ColorMap        = ColorMapB;
    ColorBrightMW   = [0.8 0.8 0.8];    % gray
    ColorDarkMW     = [1 1 1];          % white
    PlColor         = [0 0 0];          % black
    ConLinesColor   = [0 0 1];          % blue
    ConLabelsColor  = [0 0 0];          % black
    ColorIn         = [1 1 1];          % white
    ColorOut        = [1 1 1];          % white
 case 'CB'
    ColorMap        = ColorMapC;
    ColorBrightMW   = [0.0 0.0 0.15];   % blue
    ColorDarkMW     = [0 0 0];          % black
    PlColor         = [1 1 1];          % white
    ConLinesColor   = [0 1 0];          % green
    ConLabelsColor  = [0 1 1];          % cyan
    ColorIn         = [0 0 0];          % black
    ColorOut        = [0.0 0.0 0.15];   % blue
 otherwise
    %--- do nothing               ---
    %--- Use user selected colors ---
end

%---------------------------
%--- Set Planets to plot ---
%---------------------------
if (ischar(Planets)==1)
   switch Planets
    case {'no','off'}
       Planets = [];
    case 'ALL'
       Planets = (0:1:9);
    case 'All'
       Planets = [1 2 4 5 6 7 8];
    case 'EYE'
       Planets = (0:1:6);
    case 'Eye'
       Planets = [1 2 4 5 6];
    otherwise
       error('Unknown Planets option');
   end
else
   % Planets already represented as a vector of codes
end 


%---------------------------
%--- Load Star Catalogue ---
%---------------------------
switch Cat
 case 'mag6'
    load mag6.mat;
    StarCat = mag6;
    clear mag6;

    % plot nebulae only
%StarCat = NebulaeCat;


    % convert BT-VT color to B-V:
    % StarCat(:,ColCol) = 0.850.*StarCat(:,ColCol);

 case 'cats.bright.mag6'
     StarCat = cats.bright.mag6;
    
 otherwise
    load(Cat);
    
    K = findstr(Cat,'.');
    CatS = Cat(1:K-1);
    eval(['StarCat = ',CatS,';']);
end

%StarCat = [StarCat; [convertdms([03 49 14],'H','r') convertdms([+1 50 23 48],'D','R')  -2 0 -2 0 0 0]];

%---------------------------
%--- Set Magnitude Limit ---
%---------------------------
Imag    = find(StarCat(:,ColMag)<=MagLimit);
StarCat = StarCat(Imag,:);


%--------------------------------------------
%--- Applay Proper Motion to Star Catalog ---
%--------------------------------------------
switch PM
 case 'yes'
    StarCat(:,ColRA)  = StarCat(:,ColRA)  + (JD - Epoch)./JulianYear .*StarCat(:,ColPMRA)./(1000.*RAD.*ARCSEC_DEG);
    StarCat(:,ColDec) = StarCat(:,ColDec) + (JD - Epoch)./JulianYear .*StarCat(:,ColPMDec)./(1000.*RAD.*ARCSEC_DEG);
 case 'no'
    % do nothing
 otherwise
    error('Unknown PM option');
end

%----------------------------
%--- Precess Star Catalog ---
%----------------------------
switch Precess
 case 'yes'
    CatCosDir   = cosined(StarCat(:,[ColRA, ColDec]));

    RotMatP1   = rotm_coo('p',Equinox);
    RotMatP2   = rotm_coo('P',JD);

    PCatCosDir = cosined((RotMatP1*RotMatP2*CatCosDir.').');
    StarCat(:,ColRA)  = PCatCosDir(:,1);
    StarCat(:,ColDec) = PCatCosDir(:,2);
 case 'no'
    % do nothing
 otherwise
    error('Unknown Precess option');
end


%---------------------------------
%--- Load Constellations Lines ---
%---------------------------------
switch ConLines
 case 'yes'
    load('ConstellationLines.mat');
    CL = ConstellationLines;
    CL = CL.Cat;
    clear ConstellationLines;
    %--- convert coordinates to radians ---
    CL_RA1  = CL(:,1); %celestial.coo.convertdms(CL(:,1:3),  'H','r');
    CL_Dec1 = CL(:,2); %celestial.coo.convertdms(CL(:,4:7),  'D','R');
    CL_RA2  = CL(:,3); %celestial.coo.convertdms(CL(:,8:10), 'H','r');
    CL_Dec2 = CL(:,4); %celestial.coo.convertdms(CL(:,11:14),'D','R');

    %--- No proper motion handling!!! ---

    %--- Precess coordinates ---
    switch Precess
     case 'yes'
        CL_CosDir1 = cosined([CL_RA1, CL_Dec1]);
        CL_CosDir2 = cosined([CL_RA2, CL_Dec2]);

        RotMat    = rotm_coo('P',JD);
        PCL_Coo1  = cosined((RotMat*CL_CosDir1.').');
        PCL_Coo2  = cosined((RotMat*CL_CosDir2.').');
        CL_RA1    = PCL_Coo1(:,1);
        CL_Dec1   = PCL_Coo1(:,2);
        CL_RA2    = PCL_Coo2(:,1);
        CL_Dec2   = PCL_Coo2(:,2);
     case 'no'
        % do nothing
     otherwise
        error('Unknown Precess option');
    end

 case 'no'
    % do nothing
 otherwise
    error('Unknown ConLines option');
end
    

%---------------------
%--- Load MilkyWay ---
%---------------------
switch MW
 case 'yes'
    load('MilkyWay.mat');

    %--- Precess coordinates ---
    switch Precess
     case 'yes'
        RotMat    = rotm_coo('P',JD);

        %--- MW bright patches ---
        Len = length(MilkyWay.Bright);
        for Ip=1:1:Len
           MW_CosDir = cosined(MilkyWay.Bright{Ip});
           PMW_Coo   = cosined((RotMat*MW_CosDir.').');
           MilkyWay.Bright{Ip} = PMW_Coo;
        end
        %--- LMC ---
        MW_CosDir = cosined(MilkyWay.LMC);
        PMW_Coo   = cosined((RotMat*MW_CosDir.').');
        MilkyWay.LMC = PMW_Coo;
        %--- SMC ---
        MW_CosDir = cosined(MilkyWay.SMC);
        PMW_Coo   = cosined((RotMat*MW_CosDir.').');
        MilkyWay.SMC = PMW_Coo;
        %--- MW dark patches ---
        Len = length(MilkyWay.Dark);
        for Ip=1:1:Len
           MW_CosDir = cosined(MilkyWay.Dark{Ip});
           PMW_Coo   = cosined((RotMat*MW_CosDir.').');
           MilkyWay.Dark{Ip} = PMW_Coo;
        end

     case 'no'
        % do nothing
     otherwise
        error('Unknown Precess option');
    end

 case 'no'
    % do nothing
 otherwise
    error('Unknown MW option');
end



%----------------------------
%--- Plot Monthly Sky Map ---
%----------------------------
%figure(1);
set(gcf,'Visible',Visible);
%--- Plot Reference horizon circle ---
Hcirc=plot.plot_ellipse([0 0],[1 1],0,0,'k',2,1,ColorIn);
hold on;

%---------------------
%--- Plot MilkyWay ---
%---------------------
switch MW
 case 'yes'
    %--- MW bright patches ---
    Poly = MilkyWay.Bright;
    for I=1:1:length(Poly)
       MW_Cat = Poly{I};
       MW_HorizCoo = horiz_coo(MW_Cat(:,1:2),JD,GeodPos,'h');
       Ref         = refraction(MW_HorizCoo(:,2));
       MW_HorizCoo(:,2) = MW_HorizCoo(:,2) + Ref;
       [MaxAlt]    = max(MW_HorizCoo(:,2));
       if (MaxAlt<0)
          %--- do not plot segment
       else
          N_MW   = size(MW_HorizCoo,1);
          MW_HorizCat   = [MW_HorizCoo(:,1)+pi./2,pi./2-MW_HorizCoo(:,2),ones(N_MW,2)];
          [X,Y] = pr_stereographic_polar(MW_HorizCat(:,1),MW_HorizCat(:,2));
          H = patch(X,Y,'b');
          set(H,'FaceColor',ColorBrightMW,'EdgeColor',ColorBrightMW);
       end
    end
    %--- LMC ---
    Poly = MilkyWay.LMC;
    MW_Cat = Poly;
    MW_HorizCoo = horiz_coo(MW_Cat(:,1:2),JD,GeodPos,'h');
    Ref         = refraction(MW_HorizCoo(:,2));
    MW_HorizCoo(:,2) = MW_HorizCoo(:,2) + Ref;
    [MaxAlt]    = max(MW_HorizCoo(:,2));
    if (MaxAlt<0)
       %--- do not plot segment
    else
       N_MW   = size(MW_HorizCoo,1);
       MW_HorizCat   = [MW_HorizCoo(:,1)+pi./2,pi./2-MW_HorizCoo(:,2),ones(N_MW,2)];
       [X,Y] = pr_stereographic_polar(MW_HorizCat(:,1),MW_HorizCat(:,2));
       H = patch(X,Y,'b');
       set(H,'FaceColor',ColorBrightMW,'EdgeColor',ColorBrightMW);
    end
    %--- SMC ---
    Poly = MilkyWay.SMC;
    MW_Cat = Poly;
    MW_HorizCoo = horiz_coo(MW_Cat(:,1:2),JD,GeodPos,'h');
    Ref         = refraction(MW_HorizCoo(:,2));
    MW_HorizCoo(:,2) = MW_HorizCoo(:,2) + Ref;
    [MaxAlt]    = max(MW_HorizCoo(:,2));
    if (MaxAlt<0)
       %--- do not plot segment
    else
       N_MW   = size(MW_HorizCoo,1);
       MW_HorizCat   = [MW_HorizCoo(:,1)+pi./2,pi./2-MW_HorizCoo(:,2),ones(N_MW,2)];
       [X,Y] = pr_stereographic_polar(MW_HorizCat(:,1),MW_HorizCat(:,2));
       H = patch(X,Y,'b');
       set(H,'FaceColor',ColorBrightMW,'EdgeColor',ColorBrightMW);
    end	
    %--- MW dark patches ---
    Poly = MilkyWay.Dark;
    for I=1:1:length(Poly)
       MW_Cat = Poly{I};
       MW_HorizCoo      = horiz_coo(MW_Cat(:,1:2),JD,GeodPos,'h');
       Ref              = refraction(MW_HorizCoo(:,2));
       MW_HorizCoo(:,2) = MW_HorizCoo(:,2) + Ref;
       [MaxAlt]    = max(MW_HorizCoo(:,2));
       if (MaxAlt<0)
          %--- do not plot segment
       else
          N_MW   = size(MW_HorizCoo,1);
          MW_HorizCat   = [MW_HorizCoo(:,1)+pi./2, pi./2-MW_HorizCoo(:,2), ones(N_MW,2)];
          [X,Y] = pr_stereographic_polar(MW_HorizCat(:,1), MW_HorizCat(:,2));
          H = patch(X,Y,'b');
          set(H,'FaceColor',ColorDarkMW,'EdgeColor',ColorDarkMW);
       end
    end
 case 'no'
    % do nothing
 otherwise
    error('Unknown MW option');
end



%---------------------------------
%--- Plot Constellations Lines ---
%---------------------------------
switch ConLines
 case 'yes'
    CL_HorizCoo1 = horiz_coo([CL_RA1, CL_Dec1],JD,GeodPos,'h');
    CL_HorizCoo2 = horiz_coo([CL_RA2, CL_Dec2],JD,GeodPos,'h');
    Ref1         = refraction(CL_HorizCoo1(:,2));
    Ref2         = refraction(CL_HorizCoo2(:,2));
    CL_HorizCoo1(:,2) = CL_HorizCoo1(:,2) + Ref1;
    CL_HorizCoo2(:,2) = CL_HorizCoo2(:,2) + Ref2;
    I            = find(CL_HorizCoo1(:,2)>AltLimit | CL_HorizCoo2(:,2)>AltLimit);
    CL_HorizCoo1 = CL_HorizCoo1(I,:);
    CL_HorizCoo2 = CL_HorizCoo2(I,:);
    CL_HorizCat1 = [CL_HorizCoo1(:,1)+pi./2, pi./2-CL_HorizCoo1(:,2)];
    CL_HorizCat2 = [CL_HorizCoo2(:,1)+pi./2, pi./2-CL_HorizCoo2(:,2)];
    [X1,Y1] = pr_stereographic_polar(CL_HorizCat1(:,1),CL_HorizCat1(:,2));
    [X2,Y2] = pr_stereographic_polar(CL_HorizCat2(:,1),CL_HorizCat2(:,2));

    N_CL = length(I);
    for Icl=1:1:N_CL
       Hcl = plot([X1(Icl);X2(Icl)],[Y1(Icl);Y2(Icl)],'-');
       set(Hcl,'Color',ConLinesColor,'LineWidth',ConLinesWidth);
    end

 case 'no'
    % do nothing
 otherwise
    error('Unknown ConLines option');
end



%---------------------
%--- Plot Ecliptic ---
%---------------------
switch lower(Ecliptic)
 case {'yes','y','on'}
    Ec.EcLong   = [0:1:360]';
    Ec.EcLat    = zeros(length(Ec.EcLong),1);
    % convert to J2000.0 equatorial coordinates (change to equinox of date
    % in future versions).
    Ec.Eq       = coco([Ec.EcLong, Ec.EcLat]./RAD,'e','j2000.0');

    Ec.HorizCoo = horiz_coo(Ec.Eq,JD,GeodPos,'h');
    Ec.Ref      = refraction(Ec.HorizCoo(:,2));
    Ec.HorizCoo(:,2) = Ec.HorizCoo(:,2) + Ec.Ref;
    Ec.I        = find(Ec.HorizCoo(:,2)<AltLimit);
    Ec.HorizCoo(Ec.I,2) = NaN;
    Ec.HorizCat = [Ec.HorizCoo(:,1)+pi./2,pi./2-Ec.HorizCoo(:,2),Ec.Eq(:,:)];

    [Ec.X1,Ec.Y1] = pr_stereographic_polar(Ec.HorizCat(:,1),Ec.HorizCat(:,2));
    plot(Ec.X1,Ec.Y1,'k-','Color',[0.5 0.5 0.5]);

 otherwise
    % donot plot ecliptic
end


%------------------
%--- Plot Stars ---
%------------------
HorizCoo = horiz_coo(StarCat(:,[ColRA,ColDec]),JD,GeodPos,'h');
Ref      = refraction(HorizCoo(:,2));
HorizCoo(:,2) = HorizCoo(:,2) + Ref;
I        = find(HorizCoo(:,2)>AltLimit);
HorizCat = [HorizCoo(I,1)+pi./2,pi./2-HorizCoo(I,2),StarCat(I,3:4)];

switch StarType
 case 'black-edge'
    plot_smap(HorizCat,'Center',[0 0],'ColColor',4,'Project','stereographic_polar','MagFun',MagFun,'ColorMap',ColorMap,'PlotTitle','no');
 case 'same-edge'
    plot_smap(HorizCat,'Center',[0 0],'ColColor',4,'Project','stereographic_polar','MagFun',MagFun,'ColorMap',ColorMap,'PlotTitle','no','ColTypeE','0');
 otherwise
    error('Unknown StarType option');
end

H_MainAxis = gca;
set(H_MainAxis,'XDir','normal');
axis equal
AxisPosition = [0.05 0.05 0.9 0.9];
set(H_MainAxis,'Position',AxisPosition);


%--------------------
%--- Plot Planets ---
%--------------------
Nplanets = length(Planets);
Iplcoo   = 0;
PlCoo    = zeros(0,3);
for Ipl=1:1:Nplanets
   switch Planets(Ipl)
    case 0
       %--- Sun ---
       [SunRA, SunDec]  = suncoo(PlJD,'a');
       PlCoo            = [SunRA, SunDec, ones(length(SunRA),1).*-26]; 
    case 3
       %--- Moon ---
       [MoonRA,MoonDec] = mooncool(PlJD+DeltaT,GeodPos,'b');
       PlCoo            = [MoonRA, MoonDec, ones(length(MoonRA),1).*-5]; 
    case 1
       %--- Mercury ---
       %[Coo,Dist,Ang,Mag] = planet_ephem(PlJD+DeltaT,'Mercury','Topo','date',GeodPos);
       [Coo,Dist,Ang,Mag] = planet_lowephem(PlJD+DeltaT,'Mercury','Earth','SphericEq','date');
       PlCoo              = [Coo(:,1:2), Mag(:,1)];
    case 2
       %--- Venus ---
       %[Coo,Dist,Ang,Mag] = planet_ephem(PlJD+DeltaT,'Venus','Topo','date',GeodPos);
       [Coo,Dist,Ang,Mag] = planet_lowephem(PlJD+DeltaT,'Venus','Earth','SphericEq','date');
       PlCoo              = [Coo(:,1:2), Mag(:,1)];
    case 4
       %--- Mars ---
       %[Coo,Dist,Ang,Mag] = planet_ephem(PlJD+DeltaT,'Mars','Topo','date',GeodPos);
       [Coo,Dist,Ang,Mag] = planet_lowephem(PlJD+DeltaT,'Mars','Earth','SphericEq','date');
       PlCoo              = [Coo(:,1:2), Mag(:,1)];
    case 5
       %--- Jupiter ---
       %[Coo,Dist,Ang,Mag] = planet_ephem(PlJD+DeltaT,'Jupiter','Topo','date',GeodPos);
       [Coo,Dist,Ang,Mag] = planet_lowephem(PlJD+DeltaT,'Jupiter','Earth','SphericEq','date');
       PlCoo              = [Coo(:,1:2), Mag(:,1)];
    case 6
       %--- Saturn ---
       %[Coo,Dist,Ang,Mag] = planet_ephem(PlJD+DeltaT,'Saturn','Topo','date',GeodPos);
       [Coo,Dist,Ang,Mag] = planet_lowephem(PlJD+DeltaT,'Saturn','Earth','SphericEq','date');
       PlCoo              = [Coo(:,1:2), Mag(:,1)];
    case 7
       %--- Uranus ---
       %[Coo,Dist,Ang,Mag] = planet_ephem(PlJD+DeltaT,'Uranus','Topo','date',GeodPos);
       [Coo,Dist,Ang,Mag] = planet_lowephem(PlJD+DeltaT,'Uranus','Earth','SphericEq','date');
       PlCoo              = [Coo(:,1:2), Mag(:,1)];
    case 8
       %--- Neptune ---
       %[Coo,Dist,Ang,Mag] = planet_ephem(PlJD+DeltaT,'Neptune','Topo','date',GeodPos);
       [Coo,Dist,Ang,Mag] = planet_lowephem(PlJD+DeltaT,'Neptune','Earth','SphericEq','date');
       PlCoo              = [Coo(:,1:2), Mag(:,1)];
    otherwise
       error('Unkonwn planet code');
   end
   %--- convert to horizontal coordinates in the map epoch! ---
   [X,Y, RestCol] = coo2xy(PlCoo,JD,GeodPos);

   if (isempty(X)==1)
      %--- do nothing ---
   else 
      switch PlMarkSize
       case 'Mag'
          Hpl = plot(X,Y);
          set(Hpl,'LineStyle',PlLine);
          set(Hpl,'Marker',PlMarker);

          MagMarkSize = interp1(MagFun(:,1),MagFun(:,2),RestCol(:,1),'linear');
          set(Hpl,'MarkerSize',mean(MagMarkSize));

          set(Hpl,'MarkerFaceColor',PlColor,'MarkerEdgeColor',PlColor,'Color',PlColor);

       case 'Icon'
%          error('Icon option not available yet');
          set(gcf,'CurrentAxes',H_MainAxis);
          load Icon.mat
          PlanetSize = [0.1 0.1];
	  PlanetColorRange = [0 1];
          % convert X,Y coordinates to axes position coordinates
          % -1 -> 0.05 (AxisPosition...)
          % +1 -> 0.95
          AxisX = AxisPosition(1) + AxisPosition(3)./(2.*MapRadius) .* (X+MapRadius);
          AxisY = AxisPosition(1) + AxisPosition(3)./(2.*MapRadius) .* (Y+MapRadius);

          H=insert_image(X-PlanetSize(1), Y-PlanetSize(2), Icon.Planet{Planets(Ipl)}.Image,PlanetSize,PlanetColorRange,'cen',Icon.Planet{Planets(Ipl)}.Alpha);

       otherwise
          %set(Hpl,'MarkerSize',PlMarkSize);
      end

   end
end


%-----------------------------------------------
%--- Plot Additional Objects from ObjectList ---
%-----------------------------------------------
if (isempty(ObjectList)==1)
   %--- do nothing ---
else
   %--- ObjectList columns ---
   %[RA, Dec, Equinox, LineType, SymbolType, Mag, SymbolSize,
   % SymbolColor(3columns), TailPA, TailLength].
   %--- seprate different objects ---
   IndNan = find(isnan(ObjectList(:,1))==1);
   IndNan = [IndNan,size(ObjectList,1)+1];
   LastInd = 1;
   for Ido=1:1:length(IndNan)
      CurrObject = ObjectList(LastInd:IndNan(Ido)-1,:);
      %--- precess coordinates ---
      switch Precess
       case 'yes'
          ObjCosDir     = cosined(CurrObject(:,[1, 2]));
          %--- assuming the same Equinox for each object ---
          RotMatObjP1   = rotm_coo('p',CurrObject(1,3));
          RotMatObjP2   = rotm_coo('P',JD);

          PObjCosDir = cosined((RotMatObjP1*RotMatObjP2*ObjCosDir.').');
          CurrObject(:,1)  = PObjCosDir(:,1);
          CurrObject(:,2)  = PObjCosDir(:,2);
       case 'no'
          % do nothing
       otherwise
          error('Unknown Precess option');
      end
      %---
      [X,Y, RestCol] = coo2xy(CurrObject,JD,GeodPos);

      %--- plot current object ---
      Hcurobj = plot(X,Y);

      switch RestCol(1,2)
       case 0
          set(Hcurobj,'LineStyle','None');
       case 1
          set(Hcurobj,'LineStyle','-');
       case 2
          set(Hcurobj,'LineStyle','--');
       case 3
          set(Hcurobj,'LineStyle','-.');
       case 4
          set(Hcurobj,'LineStyle',':');
       otherwise
          error('Unknown line style option in ObjectList');
      end
      % 0-no symbol, 1-'o', 2-'^', 3-'s', 4-'h', 5-comet, 6-meteor shower
      switch RestCol(1,3)
       case 0
          set(Hcurobj,'Marker','None');
       case 1
          set(Hcurobj,'Marker','o','MarkerSize',RestCol(1,5));
       case 2
          set(Hcurobj,'Marker','^','MarkerSize',RestCol(1,5));
       case 3
          set(Hcurobj,'Marker','s','MarkerSize',RestCol(1,5));
       case 4
          set(Hcurobj,'Marker','p','MarkerSize',RestCol(1,5));
       case 5
          set(Hcurobj,'Marker','h','MarkerSize',RestCol(1,5));
       case 6
          set(Hcurobj,'Marker','o','MarkerSize',RestCol(1,5));
          error('Comet tail option is not available in this version');
       case 7
          set(Hcurobj,'Marker','o');
          error('Meteor shower option is not available in this version');
       otherwise
          error('Unknown marker option in ObjectList');
      end
      set(Hcurobj,'Color',RestCol(1,[6 7 8]),'MarkerFaceColor',RestCol(1,[6 7 8]),'MarkerEdgeColor',RestCol(1,[6 7 8]));

   end
end



%--------------------------------
%--- Plot Constellation Names ---
%--------------------------------
[C1,C2,C3,C4,C5,C6,C7,ConstName] = textread('Const_Names.dat','%d %d %d %d %d %d %d %s');
CN_RA  = convertdms([C1 C2 C3],'H','r');
CN_Dec = convertdms([C4 C5 C6 C7],'D','R');
[X,Y, RestCol] = coo2xy([CN_RA, CN_Dec, [1:1:length(C1)]' ],JD,GeodPos);
N_CN   = length(X);
switch ConLabels
 case 'yes'
    for Icn=1:1:N_CN
       Htext = text(X(Icn),Y(Icn),ConstName{RestCol(Icn)});
       set(Htext,'FontSize',8,'Color',ConLabelsColor)
       switch ConLabelsOrient
        case 'south'
           %--- do nothing ---
        case 'radial'
           PA_Label = atan2(Y(Icn),X(Icn));
           set(Htext,'rotation',PA_Label.*RAD+90);
        otherwise
           error('Unknown ConLabelsOrient option');
       end
    end
 case 'no'
    % do nothing
 otherwise('Unkown ConLabels option');
end


%----------------------
%--- Replot Horizon ---
%----------------------
Xcirc=get(Hcirc,'XData');
Ycirc=get(Hcirc,'YData');

In = find(Xcirc<0);
XYn = sortrows([Xcirc(In) Ycirc(In)],2);
XYn = [-1 -1; XYn; -1 1];
H=patch(XYn(:,1),XYn(:,2),'b');
set(H,'FaceColor',ColorOut,'EdgeColor',ColorOut);

Ip = find(Xcirc>0);
XYp = sortrows([Xcirc(Ip) Ycirc(Ip)],2);
XYp = [1 -1; XYp; 1 1];
H=patch(XYp(:,1),XYp(:,2),'b');
set(H,'FaceColor',ColorOut,'EdgeColor',ColorOut);

plot(Xcirc,Ycirc,'k-');
set(gcf,'Color',ColorOut);

axis off;
axis([-1 1 -1 1]);

%--------------------
%--- Plot N/E/S/W ---
%--------------------
Offset = 0.02;
Letter = 0.05;
text(0,1+Offset,'N');
text(-1-Offset-Letter,0,'E');
text(0,-1-Offset-Letter,'S');
text(1+Offset,0,'W');


%----------------------
%--- Plot Copyright ---
%----------------------
switch Copyright
 case {'None','off'}
    CopyrightText = '';
 case 'TAU'
    CopyrightText = 'Copyrights TAU AstroClub - http://astroclub.tau.ac.il';
 case 'EO'
    CopyrightText = 'Copyrights Eran Ofek - http://wise-obs.tau.ac.il/~eran/';
 otherwise
    error('Unknown Copyright option');
end
H = text(-1.2,-1.0,CopyrightText);
set(H,'rotation',90);


%-------------------
%--- Plot Legend ---
%-------------------
Legend = 'Mag';
switch Legend
 case 'Mag'
    %--- Magnitude legend ---
    LegOffset = 0.00;
    H_LegM    = axes('position',[0.91 0.05 0.06 0.26]);
    LegendCat = [0.7 0.10 5 1    
                 0.7 0.25 4 1
                 0.7 0.40 3 1
                 0.7 0.55 2 1
                 0.7 0.70 1 1
                 0.7 0.85 0 1
                 0.7 1.00 -1 1];
    LegendStarSize = interp1(MagFun(:,1),MagFun(:,2),LegendCat(:,3));

    for Ileg=1:1:size(LegendCat,1)-1
       plot(LegendCat(Ileg,1),LegendCat(Ileg,2),'o','MarkerSize',LegendStarSize(Ileg),'MarkerFaceColor',LegendStarColor,'MarkerEdgeColor',LegendStarColor);
       hold on;
       Htextleg = text(0.2,LegendCat(Ileg,2)-LegOffset,sprintf('%d',LegendCat(Ileg,3)));
    end
    box on;
    axis([0 1 0 1]);
    set(H_LegM,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');

    %--- Spectral class legend ---
    H_LegS = axes('position',[0.84 0.05 0.06 0.26]);
    LegendSpType = ['O';'B';'A';'F';'G';'K';'M']; 
    for Ileg=1:1:size(LegendCat,1)
       plot(LegendCat(Ileg,1),LegendCat(Ileg,2).*0.9,'o','MarkerSize',8,'MarkerFaceColor',ColorMapC(Ileg,1:3),'MarkerEdgeColor',ColorMapC(Ileg,1:3));
       hold on;
       Htextleg = text(0.2,LegendCat(Ileg,2).*0.9-LegOffset,sprintf('%s',LegendSpType(Ileg)));
    end
    box on;
    axis([0 1 0 1]);
    set(H_LegS,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');

    


 case {'None','off'}
    % do nothing
 otherwise
    error('Unknown Legend Option');
end


%------------------------
%--- Plot Moon Phases ---
%------------------------
switch MoonPhases
 case {'y','yes'}

    H_LegMP    = axes('position',[0.01 0.80 0.23 0.17]);
    box on;
    set(H_LegMP,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');

    DateJD  = jd2date(JD);
    MidJD   = julday([15 DateJD(2), DateJD(3)]);
    StartJD = MidJD-45;
    EndJD   = MidJD+70;
    MP      = moon_phases(StartJD,EndJD);
    Ind1    = min(find(MP(:,1)==0));
    MP      = MP(Ind1:Ind1+11,:);

    ApproxTimeZone = floor(GeodPos(1).*RAD./15)./24;   % days
    MP(:,2) = MP(:,2) + ApproxTimeZone;   % Moon phase in approximate local time
    DateMP  = jd2date(MP(:,2));

    [MP(:,1), DateMP(:,1:2)]
    Ht = text(0.05,0.90,sprintf('%02d/%02d',DateMP(1+0,1:2)));
    Ht = text(0.05,0.75,sprintf('%02d/%02d',DateMP(1+4,1:2)));
    Ht = text(0.05,0.60,sprintf('%02d/%02d',DateMP(1+8,1:2)));

    Ht = text(0.29,0.90,sprintf('%02d/%02d',DateMP(2+0,1:2)));
    Ht = text(0.29,0.75,sprintf('%02d/%02d',DateMP(2+4,1:2)));
    Ht = text(0.29,0.60,sprintf('%02d/%02d',DateMP(2+8,1:2)));

    Ht = text(0.53,0.90,sprintf('%02d/%02d',DateMP(3+0,1:2)));
    Ht = text(0.53,0.75,sprintf('%02d/%02d',DateMP(3+4,1:2)));
    Ht = text(0.53,0.60,sprintf('%02d/%02d',DateMP(3+8,1:2)));

    Ht = text(0.77,0.90,sprintf('%02d/%02d',DateMP(4+0,1:2)));
    Ht = text(0.77,0.75,sprintf('%02d/%02d',DateMP(4+4,1:2)));
    Ht = text(0.77,0.60,sprintf('%02d/%02d',DateMP(4+8,1:2)));

    MoonImSize = [0.40 0.40];
    H=insert_image(0.14,0.25,'NewMoon.jpg',MoonImSize);
    H=insert_image(0.38,0.25,'FirstQuarter.jpg',MoonImSize);
    H=insert_image(0.62,0.25,'FullMoon.jpg',MoonImSize);
    H=insert_image(0.86,0.25,'LastQuarter.jpg',MoonImSize);



 case {'n','no'}
    % do nothing
 otherwise
    error('Unknown MoonPhases option');
end




set(gcf,'InvertHardcopy','off');




%--------------------------------------------------------
%--------------------------------------------------------
%--------------------------------------------------------
function [X,Y, RestCol] = coo2xy(CooData,JD,GeodPos)
%--------------------------------------------------------
% Given [RA, Dec, additional columns]
% calculate the x,y, position in sterographic polar projection
% and return [X, Y, Rest of colums] for objects found
% above the horizon
%--------------------------------------------------------
import celestial.coo.*
import celestial.proj.*

AltLimit          = 0;
HorizCoo          = horiz_coo(CooData(:,1:2),JD,GeodPos,'h');
Ref               = refraction(HorizCoo(:,2));
HorizCoo(:,2)     = HorizCoo(:,2) + Ref;
Iabove            = find(HorizCoo(:,2)>AltLimit);
HorizCoo          = HorizCoo(Iabove,:);
if (size(CooData,2)>2)
   RestCol        = CooData(Iabove,3:end);
else
   RestCol        = [];
end
N_data            = size(HorizCoo,1);
HorizCat          = [HorizCoo(:,1)+pi./2, pi./2-HorizCoo(:,2), ones(N_data,2)];
[X,Y]             = pr_stereographic_polar(HorizCat(:,1), HorizCat(:,2));


