function prep_dss_fc(varargin)
%-----------------------------------------------------------------------------
% prep_dss_fc function                                               AstroMap
% Description: Prepare DSS finding charts, with labels, compass, and slits.
% Input  : * Arbitrary number of pairs of parameters (...,keyword,value,...):
%            Possible keywords:
%            'ListOpt1'  - List file name containing [RA,Dec,Equinox,Name],
%                          where RA and Dec are in sexagesimal format,
%                          Equinox in Years, SlitPA in deg, SlitWidth and
%                          SlitLength in arcsec.
%            'ListOpt2'  - List file names containing:
%                          [RA,Dec,Equinox,Name,SlitPA,SlitWidth,SlitLength].
%            'ListOptW'  - List file names containing: 
%                          [#,RA,Dec,Equinox,Name,...]
%            'Inv'       - Inverse colormap, {'y' | 'n'}, default is 'y'.
%            'PlotSlit'  - {'y' | 'n'}, default is 'y' (if available).
%            'PlotObj'   - {'n','bs','ro',...}. Plot circle/box around
%                          centeral coordinates, default is 'bo'.
%            'PlotSize'  - Marker size, default is 16.
%            'PlotCirc'  - {'n' | [R G B]}. Plot circle, default is 'n'.
%            'CircRad'   - Circle radius in arcsec, default is 30.
%            'PlotName'  - Plot name on finding chart {'y' | 'n'},
%                          default is 'y'.
%            'PlotCoo'   - Plot Object coordinates {'y' | 'n'},
%                          default is 'y'.
%            'PlotComp'  - Plot Compass {'y' | 'n'}, default is 'y'.
%            'PlotScale' - Scale length [arcsec], default is NaN
%                          (don't plot scale).
%            'RA'        - Vector of RA [radians], [H M S] ore sexagesimal.
%            'Dec'       - Vector of Dec [radians], [Sign D M S], or
%                          sexagesimal.
%            'FOV'       - Field of view [arcmin arcmin], default is [12 12].
%            'Equinox'   - Vector of Equinox [Years], default is [2000],
%            'Name'      - Cell array of names.
%            'SlitPA'    - Vector of Slit PA [deg], default is [],
%            'SlitWidth' - Vector of Slit width [arcsec], default is 0
%            'SlitLength'- Vector of Slit length [arcsec], default is 200.
%            'CLF'       - Clear figure after plot {'y' | 'n'}, default is 'n'.
%            'Save'      - Save finding chart using "name":
%                          {'n' | 'jpg' | 'eps'}, default is 'jpg'.
%            'Filter'    - POSS filter/epoch, see get_dss.m for options,
%                          default is '2'.
%            'Color'     - Text/Scale/Compass color, default is 'b'.
%            'SaveFITS'  - Svae FITS image {'y' | 'n'}, default is 'n';
%            'CompassPos'- Compass relative position, default is [0.1 0.1].
%            'CompassLength'-Compass relative length, default is 0.05.
%            'CooPos'    - Coordinates text relative position,
%                          default is [0.1,0.95].
%            'CooFontSize'-Coordinates text font size, default is 20.
%            'NamePos'   - Name text position, default is [0.1,0.85].
%            'NameFontSize'-Name text font size, default is 20.
%            'ScalePos'  - Scale relative position, default is [0.8,0.1].
%            'ScaleFontSize'- Scale text font sizem, default is 18.
%            'Z1_Per'    - Z1 scaling (percentile), default is 0.50.
%            'Z2_Per'    - Z2 scaling (percentile), default is 0.99.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: prep_dss_fc('ListOpt1','FC1.list','PlotObj','b+'); 
%          Name{1}='Example';
%          prep_dss_fc('RA',[10 04 00],'Dec',[1 41 0 0],'Name',Name,'PlotScale',120,'SlitPA',45,'SlitWidth',30,'Color','c','Inv','n');
% Reliable: 2
%-----------------------------------------------------------------------------



N = length(varargin);
if (0.5.*N~=floor(0.5.*N)),
   error('Illegal number of input arguments');
end

%--- Default parameters ---
ListOpt1   = [];
ListOpt2   = [];
ListOptW   = [];
Inv        = 'y';
PlotSlit   = 'y';
PlotObj    = 'bo';
PlotSize   = 16;
PlotCirc   = 'n';
CircRad    = 30;   
PlotName   = 'y';
PlotCoo    = 'y';
PlotComp   = 'y';
PlotScale  = NaN;
RA         = [];
Dec        = [];
FOV        = [12 12];
Equinox    = 2000.0;
Name       = {''};
SlitPA     = [];
SlitWidth  = 0;
SlitLength = 200;
CLF        = 'n';
Save       = 'jpg';
Filter     = '2';
Color      = 'b';
SaveFITS      = 'n';
CompassPos    = [0.1 0.1];     % compass relative position i,j
CompassLength = 0.05;          % compass relative length
CooPos        = [0.1,0.95];    % Coordinates relative position i,j
CooFontSize   = 20;            % Coordinates Font size
NamePos       = [0.1,0.85];    % Object name relative position i,j
NameFontSize  = 20;            % Name font size
ScalePos      = [0.8,0.1];     % Scale relative position i,j
ScaleFontSize = 18;            % Scale font size
Z1_Per        = 0.50;          % Z1 percentage
Z2_Per        = 0.99;          % Z2 percentage


for I=1:2:N-1,
   switch varargin{I}
    case 'ListOpt1'
       ListOpt   = varargin{I+1};

       %--- read List /Option1 ---
       [RA,Dec,Equinox,Name] = textread(ListOpt,'%s %s %f %s\n');
       RA  = convertdms(RA,'gH','H');
       Dec = convertdms(Dec,'gD','D');


    case 'ListOpt2'
       ListOpt   = varargin{I+1};

       %--- read List /Option2 ---
       [RA,Dec,Equinox,Name,SlitPA,SlitWidth,SlitLength] = textread(ListOpt,'%s %s %f %s %f %f %f\n');
       RA  = convertdms(RA,'gH','H');
       Dec = convertdms(Dec,'gD','D');

    case 'ListOptW'
       ListOpt   = varargin{I+1};

       %--- read List /Option Wise-Obs format ---
       [Number,RA,Dec,Equinox,Name] = textread(ListOpt,'%s %s %s %f %s\n');
       RA  = convertdms(RA,'gH','H');
       Dec = convertdms(Dec,'gD','D');

    case 'Inv'
       Inv       = varargin{I+1};
    case 'PlotSlit'
       PlotSlit  = varargin{I+1};
    case 'PlotObj'
       PlotObj   = varargin{I+1};
    case 'PlotCirc'
       PlotCirc  = varargin{I+1};
    case 'CircRad'
       CircRad  = varargin{I+1};
    case 'PlotSize'
       PlotSize  = varargin{I+1};
    case 'PlotName'
       PlotName  = varargin{I+1};
    case 'PlotCoo'
       PlotCoo   = varargin{I+1};
    case 'PlotComp'
       PlotComp  = varargin{I+1};
    case 'PlotScale'
       PlotScale = varargin{I+1};
    case 'RA'
       RA        = varargin{I+1};
    case 'Dec'
       Dec       = varargin{I+1};
    case 'FOV'
       FOV       = varargin{I+1};
    case 'Equinox'
       Equinox   = varargin{I+1};
    case 'Name'
       Name      = varargin{I+1};
    case 'SlitPA'
       SlitPA    = varargin{I+1};
    case 'SlitWidth'
       SlitWidth = varargin{I+1};
    case 'SlitLength'
       SlitLength= varargin{I+1};
    case 'CLF'
       CLF       = varargin{I+1};
    case 'Save'
       Save      = varargin{I+1};
    case 'Filter'
       Filter    = varargin{I+1};
    case 'Color'
       Color     = varargin{I+1};
    case 'SaveFITS'
       SaveFITS      = varargin{I+1};
    case 'CompassPos'
       CompassPos    = varargin{I+1};
    case 'CompassLength'
       CompassLength = varargin{I+1};
    case 'CooPos'
       CooPos        = varargin{I+1};
    case 'CooFontSize'
       CooFontSize   = varargin{I+1};
    case 'NamePos'
       NamePos       = varargin{I+1};
    case 'NameFontSize'
       NameFontSize  = varargin{I+1};
    case 'ScalePos'
       ScalePos      = varargin{I+1};
    case 'ScaleFontSize'
       ScaleFontSize = varargin{I+1};
    case 'Z1_Per'
       Z1_Per        = varargin{I+1};
    case 'Z2_Per'
       Z2_Per        = varargin{I+1};
    otherwise
       error('Unknown keyword Option');
   end
end      

[Links,FitsName,Images] = get_dss(RA,Dec,FOV,SaveFITS,Filter);

FitsName{:}


Nim = length(Images);
if (length(Equinox)==1),
   Equinox = Equinox.*ones(Nim,1);
end
if (length(CircRad)==1),
   CircRad = CircRad.*ones(Nim,1);
end
if (length(SlitLength)==1),
   SlitLength = SlitLength.*ones(Nim,1);
end
if (length(SlitWidth)==1),
   SlitWidth = SlitWidth.*ones(Nim,1);
end
if (length(PlotScale)==1),
   PlotScale = PlotScale.*ones(Nim,1);
end

for Iim=1:1:Nim,
   ImSize = size(Images{Iim});
   ImCL   = err_cl(mat2vec(Images{Iim}),[Z1_Per;Z2_Per]);
   Z1     = ImCL(1,1);
   Z2     = ImCL(2,2);

   Scale  = FOV(1).*60./ImSize(2);

   figure;
   imshow(Images{Iim},[Z1,Z2]);
   hold on;
   % invert (Y direction)
   set(gca,'YDir','normal');

   switch Inv
    case 'y'
       ColorMap = colormap;
       colormap(flipud(ColorMap));
    case 'n'
       % do nothing
    otherwise
       error('Unknown Inv Option');
   end

   switch PlotObj
    case 'n'
       % do nothing
    otherwise
       Hobj = plot(0.5.*ImSize(2),0.5.*ImSize(1),PlotObj);
       set(Hobj,'MarkerSize',PlotSize);
   end

   if (isstr(PlotCirc)==1),
      switch PlotCirc
       case 'n'
          % do nothing
       otherwise
          error('Unknown PlotCirc Option');
      end
   else
      Hcirc = plot_ellipse([0.5.*ImSize(2),0.5.*ImSize(1)],CircRad(Iim).*Scale,0,0,PlotCirc,1);
   end

   switch PlotSlit
    case 'n'
       % do nothing
    case 'y'
       if (isempty(SlitPA)),
          % do nothing
       else
          SlitPos   = [0.5.*ImSize(2),0.5.*ImSize(1)];
          EdgeColor = Color;
          [Hslit]   = plot_slit(SlitPos,SlitLength(Iim),SlitWidth(Iim),SlitPA(Iim),EdgeColor);
       end
    otherwise
       error('Unknown PlotSlit Option');
   end

   switch PlotName
    case 'y'
       [Hname] = text(ImSize(1).*NamePos(1),ImSize(2).*NamePos(2),sprintf('%s',Name{Iim}));
       set(Hname,'FontSize',NameFontSize,'Color',Color,'Interpreter','none');

    case 'n'
       % do nothing
    otherwise
       error('Unknown PlotCoo Option');
   end


   switch PlotCoo
    case 'y'
       StrRA  = convertdms(RA(Iim,:),'gH','SH');
       StrDec = convertdms(Dec(Iim,:),'gD','SD');
       [Hcoo] = text(ImSize(1).*CooPos(1),ImSize(2).*CooPos(2),sprintf('%s %s (%6.1f)',StrRA,StrDec,Equinox(Iim)));
       set(Hcoo,'FontSize',CooFontSize,'Color',Color);
    case 'n'
       % do nothing
    otherwise
       error('Unknown PlotCoo Option');
   end


   switch PlotComp
    case 'y'
       Rotation = 0;
       Dir      = 1;
       Pos      = ImSize.*CompassPos;
       Length   = ImSize(1).*CompassLength;
       [Hcompass,Htext] = plot_compass(Pos,Length,Color,Rotation,Dir);
       set(Htext,'Color',Color);
       set(Hcompass,'LineWidth',2);
    case 'n'
       % do nothing
    otherwise
       error('Unknown PlotCompass Option');
   end

   if (isnan(PlotScale)),
      % do nothing
   else
      Pos         = [ScalePos(1).*ImSize(1), ScalePos(2).*ImSize(2)];
      ScaleLength = PlotScale(Iim);
      UnitsName   = 'arcsec';
      Orient      = 'h';
      [Hscale,Htext]=plot_scale(Pos,Scale,ScaleLength,Color,UnitsName,Orient);
      set(Htext,'Color',Color,'FontSize',ScaleFontSize);
      set(Hscale,'LineWidth',2,'Color',Color);
   end


   switch Save
    case 'n'
       % do nothing
    case 'jpg'
       eval(sprintf('print -djpeg90 %s.jpg',Name{Iim}));
    case 'eps'
       eval(sprintf('print -depsc2 %s.eps',Name{Iim}));
    otherwise
       error('Unknown Save Option');
   end

   switch CLF
    case 'y'
       close;
    case 'n'
       % do nothing
    otherwise
       error('Unknown CLF Option');
   end


end
