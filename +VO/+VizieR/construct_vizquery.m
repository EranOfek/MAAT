function Command=construct_vizquery(varargin)
% Constrct a query string for the Vizier cdsclient command line
% Package: VO.VizieR
% Description: Constrct a query string for the Vizier cdsclient command line
% Input  : * Arbitrary number of pairs of ...,keyword,value,... arguments.
%            The following keywords are supported:
%            'CatName' - The cdclient program for the specific catalog.
%                        See Catalog.VizieR.cdsclient_prog_names.
%                        Default is 'finducac4'.
%            'RA'      - J2000.0 RA [sexagesimal, radians or deg].
%            'Dec'     - J2000.0 Dec [sexagesimal, radians or deg].
%            'CooUnits'- Coordinate units. Default is 'deg'.
%            'Radius'  - Search radius. Default is 10.
%                        For box search this can be [X,Y] box size.
%            'RadiusUnits' - Search radius units. Default is 'arcmin'.
%                        See convert.angular for options.
%            'ObjName' - Object name. If provided, then will search by
%                        object name. Default is empty.
%            'RegionType' - Search region type: 'circ'|'box'.
%                        Default is 'circ'.
%            'MaxRecord' - Max. number of records. Default is 100000.
%            'PathPar'   - Additional arguments to pass to
%                          Catalog.VizieR.cdsclient_path.
%                          Default is {}.
%            'IncludePath' - Include path in command string.
%                          Default is true.
% Output : - Command string
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Path=VO.VizieR.construct_vizquery
%          Path=VO.VizieR.construct_vizquery('RA',90.0,'Dec',-15.2)
% Reliable: 2


DefV.CatName             = 'finducac4'; %'findwise';
DefV.RA                  = 180;
DefV.Dec                 = 0;
DefV.CooUnits            = 'deg';
DefV.Radius              = 10;   % or [x,y]
DefV.RadiusUnits         = 'arcmin';
DefV.ObjName             = '';
DefV.RegionType          = 'circ';  % 'circ' | 'box'
DefV.MaxRecord           = 100000;
DefV.PathPar             = {};
DefV.IncludePath         = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% convert J2000.0 RA/Dec to deg
if (ischar(InPar.RA))
    RA = celestial.coo.convertdms(InPar.RA,'SH','d');
else
    RA = convert.angular(InPar.CooUnits,'deg',InPar.RA);
end
if (ischar(InPar.Dec))
    Dec = celestial.coo.convertdms(InPar.Dec,'SD','D');
else
    Dec = convert.angular(InPar.CooUnits,'deg',InPar.Dec);
end

% convert search radius to arcmin
Radius = convert.angular(InPar.RadiusUnits,'arcmin',InPar.Radius);

if (InPar.IncludePath)
    Path = VO.VizieR.cdsclient_path(InPar.PathPar{:});
else
    Path = '';
end

% search by coordinate or object name
if (isempty(InPar.ObjName))
    % search by coordinates
    Command = sprintf('%s%s -c %11.6f %+11.6f',Path,InPar.CatName,RA,Dec);
else
    % search by object identifier/name
    Command = sprintf('%s%s -i %s',Path,InPar.CatName,InPar.ObjName);
end

% search region type
switch lower(InPar.RegionType)
    case 'circ'
        Command = sprintf('%s -r %f',Command,Radius);
    case 'box'
        if (numel(Radius)>1)
            Command = sprintf('%s -b %f,%f',Command,Radius(1),Radius(2));
        else
            Command = sprintf('%s -b %f',Command,Radius);
        end
    otherwise
        error('Unknown RegionType option');
end

% max records
if (~isempty(InPar.MaxRecord))
    Command = sprintf('%s -m %d',Command,InPar.MaxRecord);
end