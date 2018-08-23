function Out=footprint(Sim,varargin)
% Return the footprint and center for a list of HEAD objects.
% Package: @WorldCooSys
% Description: Return the footprint (verteces of images footprint) and
%              center for a list of HEAD objects.
% Input  : - An WCS or HEAD object or any class containing WCS
%            (e.g., SIM image).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'CCDSEC' - The CCD section. This can be an header keyword
%                       string from which to obtain the CCDSEC,
%                       or a CCDSEC vector [xmin xmax ymin ymax].
%                       Alternativelly, if empty will use the actual image
%                       size to estimate CCDSEC.
%                       Default is empty.
%            'OutUnits' - Output coordinates units {'rad'|'deg'}.
%                       Default is 'rad'.
% Output : - Structure containing the footprints long/lat.
%            The following fields are available:
%            .FootLong - Vector of footprints longitude.
%            .FootLat  - Vector of footprints latitude.
%            .CenLong  - Center longitude.
%            .CenLat   - Center latitude.
%            .GcenLong - Geometric center longitude.
%            .GcenLat  - Geometric center latitude.
%            .Radius   - Distance between geometric center and furtest
%                        corner.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Out,Sim]=footprint(S);
% Reliable: 2
%--------------------------------------------------------------------------
import celestial.coo.*

RAD = 180./pi;

%ImageField     = 'Im';
%HeaderField    = 'Header';
%FileField      = 'ImageFileName';
%MaskField      = 'Mask';
%BackImField    = 'BackIm';
%ErrImField     = 'ErrIm';
%WeightImField  = 'WeightIm';


DefV.CCDSEC           = [];   % [], or [xmin xmax ymin ymax] or string
DefV.OutUnits         = 'rad';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


switch lower(InPar.OutUnits)
    case 'rad'
        ConvUnits = 1;
    case 'deg'
        ConvUnits = RAD;
    otherwise
        error('Unknown OutUnits option');
end


% read images
%Sim = images2sim(Sim,varargin{:});
Nim = numel(Sim);

if (HEAD.ishead(Sim))
    CCDSEC = ccdsec(Sim,InPar.CCDSEC);
else
    if (isnumeric(InPar.CCDSEC) && ~isempty(InPar.CCDSEC))
        CCDSEC = InPar.CCDSEC;
    else
        error('For input WCS object CCDSEC must be non-empty numeric');
    end
end

Out = struct('FootLong',cell(Nim,1),'FootLat',cell(Nim,1),...
             'CenLong',cell(Nim,1),'CenLat',cell(Nim,1),...
             'GcenLong',cell(Nim,1),'GcenLat',cell(Nim,1));
for Iim=1:1:Nim
    % calculate footprints
    X  = CCDSEC(Iim,[1 2 2 1]);
    Y  = CCDSEC(Iim,[3 3 4 4]);
    CX = 0.5.*(CCDSEC(Iim,1)+CCDSEC(Iim,2));
    CY = 0.5.*(CCDSEC(Iim,3)+CCDSEC(Iim,4));
    %[Long,Lat] = xy2sky(Sim(Iim),[X CX].',[Y CY].');
    [Long,Lat] = xy2coo(Sim(Iim),[X CX].',[Y CY].');
    [CD1,CD2,CD3] = coo2cosined(Long,Lat);
    MeanCD1 = mean(CD1(1:4));
    MeanCD2 = mean(CD2(1:4));
    MeanCD3 = mean(CD3(1:4));
    
    [MeanLong,MeanLat] = cosined2coo(MeanCD1,MeanCD2,MeanCD3);
    
    Out(Iim).FootLong = Long(1:4).*ConvUnits;
    Out(Iim).FootLat  = Lat(1:4).*ConvUnits;
    Out(Iim).CenLong  = Long(5).*ConvUnits;
    Out(Iim).CenLat   = Lat(5).*ConvUnits;
    Out(Iim).GcenLong = MeanLong.*ConvUnits;
    Out(Iim).GcenLat  = MeanLat.*ConvUnits;
    % radius of circle (center to corner)
    Out(Iim).Radius   = max(sphere_dist_fast(MeanLong,MeanLat,Long(1:4),Lat(1:4))).*ConvUnits;
end

            

