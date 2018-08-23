function [OutPx,OutPy,HColor]=amapproj(Coo,Proj,AreaVecType,Symbol,Grid,RotationMatrix,SGrid)
% Old map projection function (not supported)
% Package: celestial.map
% Description: Plot a list of coordinates (as dots/lines) on a map with
%              a chosen projection.
% Input  : - matrix of coordinates.
%            [Longitude (rad), Latitude (rad), Magnitude, Spectral-type].
%            Magnitude and Spectral-type are optional.
%            In case that the "all sky map" option is selected, the
%            coordinates should be in range [-pi pi -pi/2 pi/2].
%            Magnitude cotrol the symbol size.
%            The spectral-type control the symbol color and it can be
%            one of the following:
%            1  -  Blueish [0 0.5 0.9]
%            2  -  Cyan
%            3  -  Green
%            4  -  Yellow
%            5  -  Orange [1 0.65 0]
%            6  -  Reddish    [1 0 0.5]
%            7  -  Black    [0 0 0]
%            8  -  Gray    [0.8 0.8 0.8]
%            9  -  Blcak 
%            * In case of seven columns matrix, the columns
%              are assumed to be [H M S DecSign D M S]
%          - projection type:
%            'a' - Aitoff (default)
%            'm' - Mollweide equal-area
%            'h' - Hammer
%            'p' - Parabolic
%            's' - Sinusoidal
%            'l' - Lambert equal-area cylindrical
%            'b' - Behrmann equal-area cylindrical
%            't' - Tristan Edwards cylindrical
%            'P' - Peters cylindrical
%            'G' - Gall Orthographic cylindrical
%            'B' - Balthasart cylindrical
%            'C' - Cassini
%            'x' - XY, no projection
%            'r' - polar
%            'g' - Gnomonic
%            'o' - Bonne
%            'S' - Stereographic
%            #.# - (number), conic projection with height #.#
%          - area vector
%            'a' - for all sky globe. (default) or
%            [Long_Min, Long_Max, Lat_Min, Lat_Max]
%          - object symbol and color
%            '.', 'o', '>', '<', 'd', '*', ...
%            default is '.'
%          - grid symbol and color
%            'N' - no grid.
%            default is 'r:'
%          - 3x3 rotation matrix of secondary grid lines.
%            default is no secondary grid. recomended only
%            with AreaVec='a'.
%          - Secondary grid format, default is 'r.'.
%            default is no secondary grid.
% Output : - The projected X position of the data points.
%          - The projected Y position of the data points.
% Tested : Matlab 5.2       
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: NOT WORKING ANYMORE
%--------------------------------------------------------------------------
RADIAN   = 180./pi;
TPI      = 2.*pi;
ColorMap = colormap;

% boundary settings
BSymbol              = '-';   % Boundary symbol type
BSymbolWidth         = 2.0;   % Boundary symbol width
% Grid settings
GSymbolWidth         = 0.5;   % Grid symbol width
GsSymbolWidth        = 2;     % secondary grid symbol
NLatGrid             = 13;     % number of latitude grid lines
NLongGrid            = 13;     % number of longitude grid lines
GLabelFontSize       = 6;
PlotLongitudePos     = 'm';  %plot long on lat=0 {m-middle,b-bot,t-top}
% Data settings
DataMarkerSize       = 20;   % data marker size
DataMarkerRange      = [20 50];   % range of data marker size (in case of magnitude)
% contour settings
ContourSymbol        ='-';
ContourSymbolWidth   = 1;
ContourLabelFontSize = 6;
ContourLabelType     = 'h';  % 'd'- deg.   'h'-hours/deg.

% magnitude labeling
PlotLabel            = 'n';
MagLabelFontSize     = 6;

HColor = zeros(6,1);

% setting the grid labels position parameters
if (Proj=='a' | Proj=='h'),
   IncX                = 0.1; % 0.2;
   IncY                = 0.05; % 0.1;
   IncF                = 1.05;
   IncXY               = 0.05;
elseif (Proj=='m' | Proj=='s' | Proj=='p'),
   IncX                = 0.20; % 0.2;
   IncY                = 0.05; % 0.1;
   IncF                = 1.03;
   IncXY               = 0.05;
elseif (Proj=='l' | Proj=='b' | Proj=='t' | Proj=='P' | Proj=='G' | Proj=='B' | Proj=='c'),
   IncX                = 0.15; % 0.2;
   IncY                = 0.05; % 0.1;
   IncF                = 1.00;
   IncXY               = 0.05;
elseif (Proj=='x'),
   IncX                = 0.02; % 0.2;
   IncY                = 0.02; % 0.1;
   IncF                = 1.00;
   IncXY               = 0.05;
else
   IncX                = 0.1; % 0.2;
   IncY                = 0.08; % 0.1;
   IncF                = 1.00;
   IncXY               = 0.01;
end
   
   

PlotMesh = 0;
if (nargin==1),
   % set default value
   Proj        = 'a';   % Aitoff projection
   AreaVecType = 'a';
   Symbol      = '.';
   Grid        = 'r-';
   % no rotation matrix (default)
   SecGrid     = 'N';
elseif (nargin==2),
   AreaVecType = 'a';
   Symbol      = '.';
   Grid        = 'r-';
   % no rotation matrix (default)
   SecGrid     = 'N';
elseif (nargin==3),
   Symbol  = '.';
   Grid    = 'r-';
   % no rotation matrix (default)
   SecGrid = 'N'; 
elseif (nargin==4),
   Grid    = 'r-';
   % no rotation matrix (default)
   SecGrid = 'N'; 
elseif (nargin==5),
   % no rotation matrix (default)
   SecGrid = 'N'; 
elseif (nargin==6),
   SecGrid = 'Y';
   GsSymbol = 'r.';
elseif (nargin==7),
   % no defaults   
elseif (nargin==8),
   % mesh ploting 
   PlotMesh = 1;
else
   error ('illigal number of arguments');
end
GSymbol = Grid;




axes;
axis('off');
hold on;

%----------------------------------
% check if data points are in range
%----------------------------------
NewCooTmp = zeros(size(Coo));

if (AreaVecType(1)=='a'),
   AreaVec = [-pi pi -pi./2 pi./2];
   Scale   = 1;
else
   AreaVec = AreaVecType;
   Scale   = 2.*pi./(AreaVec(2)-AreaVec(1));
   PlotLongitudePos = 'b';  % put longitude labels on bot.
   IncY    = 3.30;
   IncXY    = -0.08;
end

J = 0;
for I=1:1:length(Coo(:,1)),
   if (Coo(I,1)<AreaVec(2) & Coo(I,1)>AreaVec(1)),
      if (Coo(I,2)<AreaVec(4) & Coo(I,2)>AreaVec(3)),
         J = J + 1;
         NewCooTmp(J,:) = Coo(I,:);
      end
   end
end
NewCoo = zeros(J,2);
NewCoo = NewCooTmp(1:J,:);
Coo    = NewCoo;


Ncol = length(Coo(1,:));
% convert coordinate type to radians in range [-pi pi -pi/2 pi/2].
MagPlot  = 'n';  % magnitude plotting
SpecPlot = 'n';  % spectral type plotting
if (Ncol==2),
   X = Coo(:,1);
   Y = Coo(:,2);
elseif (Ncol==3),
   % magnitude case
   X = Coo(:,1);
   Y = Coo(:,2);
   Magnitude = Coo(:,3);
   MagPlot = 'y';
elseif (Ncol==4),
   % magnitude case
   X = Coo(:,1);
   Y = Coo(:,2);
   Magnitude = Coo(:,3);
   SpecType  = Coo(:,4);
   MagPlot  = 'y';
   SpecPlot = 'y';
elseif (Ncol==7),
   % convert HMS/DMS to deg.
   X = Coo(:,1).*15 + 15.*Coo(:,2)./60 + 15.*Coo(:,3)./3600;
   Y = Coo(:,4).*(Coo(:,5) + Coo(:,6)./60 + Coo(:,7)./3600);
   % convert to radians
   X = X./RADIAN;
   Y = Y./RADIAN;
else
   error('Illigal coordinate type');
end


%------------------------------------------------------------
% all coordinate are in radians in range [-pi pi -pi/2 pi/2].
%------------------------------------------------------------

%-------------------------------
% prepare & plot boundary vector
%-------------------------------
if (AreaVecType(1)=='a'),
   % plot all sky globe in the wanted projection
   AreaVec   = [-pi pi -pi./2 pi./2];
   PlotAllSky = 'y';
else
   % plot limitted sky map
   AreaVec = AreaVecType;
   PlotAllSky = 'n';
end
Nelements = 40;    % half number of elements in boundry vector
Nstep     = Nelements - 1;
SpanX = AreaVec(2) - AreaVec(1);
SpanY = AreaVec(4) - AreaVec(3);

IncX = IncX.*SpanX./TPI;
IncY = IncY.*SpanY./pi;

BVecX = zeros(Nelements,1);
BVecY = zeros(Nelements,1);

if (ischar(Proj) & Proj~='r'),
   % if not conic section, plot eqal longitude lines.
   BVecX = AreaVec(1).*ones(Nelements,1);
   BVecY = [AreaVec(3):(SpanY./Nstep):AreaVec(4)]';
   [Px,Py]=projectcoo(BVecX,BVecY,Proj,Scale);
   H=plot(Px,Py,BSymbol);
   set(H,'LineWidth',BSymbolWidth,'MarkerSize',BSymbolWidth);

   BVecX = AreaVec(2).*ones(Nelements,1);
   BVecY = [AreaVec(3):(SpanY./Nstep):AreaVec(4)]';
   [Px,Py]=projectcoo(BVecX,BVecY,Proj,Scale);
   H=plot(Px,Py,BSymbol);
   set(H,'LineWidth',BSymbolWidth,'MarkerSize',BSymbolWidth);
end

BVecX = [AreaVec(1):(SpanX./Nstep):AreaVec(2)]';
BVecY = AreaVec(3).*ones(Nelements,1);
[Px,Py]=projectcoo(BVecX,BVecY,Proj,Scale);
H=plot(Px,Py,BSymbol);
set(H,'LineWidth',BSymbolWidth,'MarkerSize',BSymbolWidth);

BVecX = [AreaVec(1):(SpanX./Nstep):AreaVec(2)]';
BVecY = AreaVec(4).*ones(Nelements,1);
[Px,Py]=projectcoo(BVecX,BVecY,Proj,Scale);
H=plot(Px,Py,BSymbol);
set(H,'LineWidth',BSymbolWidth,'MarkerSize',BSymbolWidth);


%-----------------------------------------------------------
% plot the data points/contour/mesh/surf/countour on sky map
%-----------------------------------------------------------
if (PlotMesh==0),
   %-----------------------
   % plot data point on map
   %-----------------------
   [Px,Py]=projectcoo(X,Y,Proj,Scale);
   if (nargout>0),
      OutPx = Px;
      OutPy = Py;
   end
   if (MagPlot=='y'),
      MinMag = min(Magnitude);
      MaxMag = max(Magnitude);
      if ((MaxMag-MinMag)~=0),
         DataMarkerVec = DataMarkerRange(1)+(DataMarkerRange(2)-DataMarkerRange(1)).*(MaxMag - Magnitude)./(MaxMag - MinMag);
      else
         DataMarkerVec = mean(DataMarkerRange).*ones(size(Magnitude));
      end
      % magnitude plotting
      % handle vector
      Hd = zeros(size(Magnitude));
      if (SpecPlot=='n'),
         for I=1:1:length(Magnitude),
            Hd(I) = plot(Px(I),Py(I),Symbol);
            set(Hd(I),'MarkerSize',DataMarkerVec(I));
            if (PlotLabel=='y'),
              if (rand(1)>0.0),
	       Ht = text(Px(I)+0.002.*Scale ,Py(I)+0.002.*Scale,num2str(Magnitude(I)));
               set(Ht,'FontSize',MagLabelFontSize);
              end
            end
         end
      else
         % plot spectral-type (by color) and magnitude (by size)
         for I=1:1:length(Magnitude),
            Hd(I) = plot(Px(I),Py(I),Symbol);
            switch SpecType(I)
               case (1),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color',[0 0.5 0.9]); % blueish
                  HColor(1) = Hd(I);
               case (2),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color','Cyan');
                  HColor(2) = Hd(I);
               case (3),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color','Green');
                  HColor(3) = Hd(I);
               case (4),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color','Yellow');
                  HColor(4) = Hd(I);
               case (5),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color',[1 0.65 0]); % orange
                  HColor(5) = Hd(I);
               case (6),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color',[1 0 0.5]); % reddish
                  HColor(6) = Hd(I);

               case (7),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color',[0 0 0]); % black
                  HColor(7) = Hd(I);

               case (8),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color',[0.8 0.8 0.8]); % gray
                  HColor(3) = Hd(I);

               case (9),
                  set(Hd(I),'MarkerSize',DataMarkerVec(I),'Color',[0 0 0],'Marker','^','MarkerSize',5); % black ^
                  HColor(7) = Hd(I);

               otherwise
                  error('Unknown color in spectral-type');
            end
         end         
      end
   else
      % no magnitude plotting
      H = plot(Px,Py,Symbol);
      set(H,'MarkerSize',DataMarkerSize);
   end
else
   %--------------
   % plot mesh map
   % (PlotMesh==1)
   %--------------
   %error('Not Yet Implemented');   
   
   % contour plot
   LongGrid = [-pi:5./RADIAN:pi];
   LatGrid  = [-pi./2:5./RADIAN:pi./2]';
   ValMat   = LatGrid*ones(1,length(LongGrid));
   ConVal   = 10;
   if (length(ConVal)==1),
      ConNum = ConVal;
   else
      ConNum = length(ConVal);
   end
   ConMat   = contourc(LongGrid, LatGrid, ValMat, ConVal);
   ConLevel = zeros(ConNum,1);
   ContColor = zeros(ConNum,3);
   CInd = 1;
   for LInd=1:1:ConNum,
      ConLevel(LInd) = ConMat(1,CInd);
      PairNum        = ConMat(2,CInd);
      Xc             = [ConMat(1,CInd+1:CInd+PairNum)]';
      Yc             = [ConMat(2,CInd+1:CInd+PairNum)]';
      [Px,Py]         = projectcoo(Xc,Yc,Proj,Scale);
      Hc(LInd)        = plot(Px,Py,ContourSymbol);
      ColorMapLine    = floor(LInd.*floor(length(ColorMap(:,1))./ConNum)) + 1;
      ContColor(LInd,:) = ColorMap(ColorMapLine,1:3);
      set(Hc(LInd),'LineWidth',ContourSymbolWidth,'Color',ContColor(LInd,:));
 
      NewConMat = [ConLevel(LInd), [Px']; PairNum, [Py']];
      Hcl  = clabel(NewConMat,Hc(LInd));
      set(Hcl,'FontSize',ContourLabelFontSize);
      CInd = CInd + PairNum + 1; 
   end
   
end


%---------------------------------------
% plot grid lines, using user parameters
%---------------------------------------
if (Grid~='N'),
   %----------------
   % plot grid lines
   %----------------
   GStepLat  = NLatGrid - 1;
   GStepLong = NLongGrid - 1;

   GVecX = zeros(Nelements,1);
   GVecY = zeros(Nelements,1);

   % plot equal latitude grid lines
   Step = SpanY./GStepLat;
   for JLat = AreaVec(3)+Step:Step:AreaVec(4)-Step,
      GVecX = [AreaVec(1):(SpanX./Nstep):AreaVec(2)]';
      GVecY = JLat.*ones(size(GVecX));
      [Px,Py]=projectcoo(GVecX,GVecY,Proj,Scale);
      H=plot(Px,Py,GSymbol);
      set(H,'LineWidth',GSymbolWidth,'MarkerSize',GSymbolWidth);
      % plot grid Latitude label
      if (~ischar(Proj) | Proj=='r'),
         % conic plots
         text(max(Px)+IncX,IncY,sprintf('%4.1f',JLat.*RADIAN),'FontSize',GLabelFontSize);
      else
         text(max(Px)+IncX,max(abs(Py)).*IncF.*sign(JLat),sprintf('%4.1f',JLat.*RADIAN),'FontSize',GLabelFontSize);
      end
   end
   % plot equal longitude grid lines
   Step = SpanX./GStepLong;
   for JLong = AreaVec(1)+Step:Step:AreaVec(2)-Step,
      GVecY = [AreaVec(3):(SpanY./Nstep):AreaVec(4)]';
      GVecX = JLong.*ones(size(GVecY));
      [Px,Py]=projectcoo(GVecX,GVecY,Proj,Scale);
      H=plot(Px,Py,GSymbol);
      set(H,'LineWidth',GSymbolWidth,'MarkerSize',GSymbolWidth);
      % plot grid Longitude label
      if (~ischar(Proj) | Proj=='r'),
         % conic plots
         % to be implemented
      else
         if (PlotLongitudePos=='m'),
            % plot longitude labels in lat=0
            if (ContourLabelType=='d'),
               text(sign(max(Px)).*max(abs(Px))+IncXY,IncY,sprintf('%4.1f',JLong.*RADIAN),'FontSize',GLabelFontSize);  
            elseif (ContourLabelType=='h'),
               if (JLong<0),
                  HourLong = (2.*pi+JLong).*RADIAN./15;
               else
                  HourLong = JLong.*RADIAN./15;
               end
               text(sign(max(Px)).*max(abs(Px))+IncXY,IncY,sprintf('%4.1f',HourLong),'FontSize',GLabelFontSize);  
            else
               error('Illegal ContourLabelType');
            end
         elseif (PlotLongitudePos=='b'),
            % plot longitude label on botoum
            if (ContourLabelType=='d'),
               text(sign(max(Px)).*max(abs(Px))+IncXY,min(Py)-IncY,sprintf('%4.1f',JLong.*RADIAN),'FontSize',GLabelFontSize);      
            elseif (ContourLabelType=='h'),
               if (JLong<0),
                  HourLong = (2.*pi+JLong).*RADIAN./15;
               else
                  HourLong = JLong.*RADIAN./15;
               end
               text(sign(max(Px)).*max(abs(Px))+IncXY,min(Py)-IncY,sprintf('%4.1f',HourLong),'FontSize',GLabelFontSize);      
            else
               error('Illegal ContourLabelType');
            end
         elseif (PlotLongitudePos=='t'),
            % plot longitude label on top
            if (ContourLabelType=='d'),
               text(sign(max(Px)).*max(abs(Px))+IncXY,max(Py)+IncY,sprintf('%4.1f',JLong.*RADIAN),'FontSize',GLabelFontSize);             
            elseif (ContourLabelType=='h'),
               if (JLong<0),
                  HourLong = (2.*pi+JLong).*RADIAN./15;
               else
                  HourLong = JLong.*RADIAN./15;
               end
               text(sign(max(Px)).*max(abs(Px))+IncXY,max(Py)+IncY,sprintf('%4.1f',HourLong),'FontSize',GLabelFontSize);             
            else
               error('Illegal ContourLabelType');
            end
         else
            error('unknown plot longitude position option, {m,b,t}');
         end
      end      
   end
else
   % user didn't ask for main grid lines
end


if (SecGrid~='N'),
   %-----------------------------------------------------
   % plot secondary grid lines, rotated to the main lines
   %-----------------------------------------------------
   GStepLat  = NLatGrid - 1;
   GStepLong = NLongGrid - 1;
   
   Nstep = 100;
   GVecX = zeros(Nelements,1);
   GVecY = zeros(Nelements,1);

   % plot equal latitude grid lines
   Step = SpanY./GStepLat;
   for JLat = AreaVec(3):Step:AreaVec(4),
      GVecX = [AreaVec(1):(SpanX./Nstep):AreaVec(2)]';
      GVecY = JLat.*ones(size(GVecX));
      % convert Long & Lat to cosine direction
      CosineDir = cosined([GVecX,GVecY]);
      % rotate cosine direction vector
      RotCosineDir = [RotationMatrix*CosineDir']';
      % convert back to long & lat
      LongLat = cosined(RotCosineDir);
      GVecX = LongLat(:,1);
      GVecY = LongLat(:,2);
      
      [Px,Py]=projectcoo(GVecX,GVecY,Proj,Scale);
      H=plot(Px,Py,GsSymbol);
      set(H,'LineWidth',GsSymbolWidth,'MarkerSize',GsSymbolWidth);
   end

   % plot equal longitude grid lines
   Step = SpanX./GStepLong;
   for JLong = AreaVec(1):Step:AreaVec(2),
      GVecY = [AreaVec(3):(SpanY./Nstep):AreaVec(4)]';
      GVecX = JLong.*ones(size(GVecY));
      % convert Long & Lat to cosine direction
      CosineDir = cosined([GVecX,GVecY]);
      % rotate cosine direction vector
      RotCosineDir = [RotationMatrix*CosineDir']';
      % convert back to long & lat
      LongLat = cosined(RotCosineDir);
      GVecX = LongLat(:,1);
      GVecY = LongLat(:,2);

      [Px,Py]=projectcoo(GVecX,GVecY,Proj,Scale);
      H=plot(Px,Py,GsSymbol);
      set(H,'LineWidth',GsSymbolWidth,'MarkerSize',GsSymbolWidth);
   end   
else
   % user didn't ask for secondary grid lines
end
   
% set the Aspect Ratio
Hax = get(gcf,'Children');
AspectRatioVec = [1 1 1];
set(Hax,'DataAspectRatio',AspectRatioVec);

% end of program
hold off;
