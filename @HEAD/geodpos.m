function Out=geodpos(Sim,varargin)
%--------------------------------------------------------------------------
% geodpos function                                             class/@HEAD
% Description: Get observatory geodetic position from image header.
% Input  : - An HEAD object (or e.g., a SIM object).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Long'     - Longitide or a string containing
%                         the longitude header keyword name.
%                         Default is 'OBSLON'.
%            'Lat'      - Latitude or a string containing
%                         the latitude header keyword name.
%                         Default is 'OBSLAT'.
%            'Height'   - Height or a string containing
%                         the height header keyword name.
%                         Default is NaN.
%            'RefEllip' - String containing the reference ellipsoid name.
%                         Default is 'WGS84'.
%            'InAngUnits'- Units of the input long/lat.
%                         Default is 'deg'.
%            'InHeightUnits' - Units of of the input height.
%                         Default is 'm'.
% Output : - A structure arry containing the Geodetic position as obtained
%            from the header. The following fields are available:
%            .RefEllip
%            .Height
%            .Long     [rad]
%            .Lat      [rad]
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Out=geodpos(S)
% Reliable: 2
%--------------------------------------------------------------------------


RAD    = 180./pi;  % 1 radian = 57.295... deg
InvRAD = pi./180;

% SIM class fields
% ImageField      = 'Im';
% HeaderField     = 'Header';
% BackImField     = 'BackIm';
% ErrImField      = 'ErrIm';
% CatField        = 'Cat';
% CatColField     = 'Col';
% CatColCellField = 'ColCell';


DefV.Long              = 'OBSLON';
DefV.Lat               = 'OBSLAT';
DefV.Height            = NaN;
DefV.RefEllip          = 'WGS84';
%DefV.ReadImage         = false;
DefV.InAngUnits        = 'deg';   % deg|rad
DefV.InHeightUnits     = 'm';     % cm|m|km
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

AngConv    = convert.angular(DefV.InAngUnits,'rad');
HeightConv = convert.units(DefV.InHeightUnits,'cm');

Out.RefEllip = InPar.RefEllip;

Keys2get = cell(0,1);
KeysName = cell(0,1);

if (isnumeric(InPar.Long)),
    Out.Long    = InPar.Long;
else
    KeysName{end+1} = 'Long';
    Keys2get{end+1} = InPar.Long;
end
if (isnumeric(InPar.Lat))
    Out.Lat     = InPar.Lat;
else
    KeysName{end+1} = 'Lat';
    Keys2get{end+1} = InPar.Lat;    
end
if (isnumeric(InPar.Height))
    Out.Height  = InPar.Height;
else
    KeysName{end+1} = 'Height';
    Keys2get{end+1} = InPar.Height;    
end

Nim = numel(Sim); % Na'ama, 20180808

if (~isempty(Keys2get))
    % read information from header
    %Sim = images2sim(Sim,varargin{:}); %,'ReadImage',InPar.ReadImage);
    %Nim = numel(Sim); % commented by Na'ama, 20180808
    
    if (~any(strcmp(KeysName,'Long')))
        [Out(1:1:Nim).Long] = deal(InPar.Long);
    end
    if (~any(strcmp(KeysName,'Lat')))
        [Out(1:1:Nim).Lat] = deal(InPar.Lat);
    end
    if (~any(strcmp(KeysName,'Height')))
        [Out(1:1:Nim).Height] = deal(InPar.Height);
    end
    
    %[~,Struct] = sim_getkeyvals(Sim,Keys2get,varargin{:});
    [CellGeodPos] = mgetkey(Sim,Keys2get);
    Struct        = cell2struct(CellGeodPos,Keys2get,2);
    for Ik=1:1:numel(KeysName),
        [Out(1:1:Nim).(KeysName{Ik})] = Struct.(Keys2get{Ik});
    end
    
    
end

for Iim=1:1:Nim,
    Out(Iim).Long   = Out(Iim).Long.*AngConv;
    Out(Iim).Lat    = Out(Iim).Lat.*AngConv;
    Out(Iim).Height = Out(Iim).Height.*HeightConv;
end 