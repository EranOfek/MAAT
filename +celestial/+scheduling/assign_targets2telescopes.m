function ResAs=assign_targets2telescopes(NewRA,NewDec,CurRA,CurDec,varargin)
%




RAD = 180./pi;

InPar = inputParser;
addOptional(InPar,'Lon',35./RAD);
addOptional(InPar,'Lat',30./RAD);
addOptional(InPar,'Ntel',8);
addOptional(InPar,'AllTelAzAltLimit',{[0 15; 90 15; 180 15; 270 15; 360 15]});
addOptional(InPar,'AllTelRangeHA',{[-120 120]});
addOptional(InPar,'AllTelRangeDec',{[-30 90]});

addOptional(InPar,'LST',[]);  % [rad]
addOptional(InPar,'InterpMethod','linear'); 

parse(InPar,varargin{:});
InPar = InPar.Results;

if isempty(InPar.LST)
    % get current LST
    JD  = celestial.time.julday;
    InPar.LST = celestial.time.lst(JD,InPar.Lon,'m');
end

NewHA = InPar.LST - NewRA;

Nazalt = numel(InPar.AllTelAzAltLimit);
Nha    = numel(InPar.AllTelRangeHA);
Ndec   = numel(InPar.AllTelRangeDec);

% check which telescope is allowed to observe each one of the new targets
Ntarget = numel(NewRA);
FlagValid = false(Ntarget,InPar.Ntel);
for Itel=1:1:InPar.Ntel
    Itel_azalt = min(Nazalt,Itel);
    Itel_ha    = min(Nha,Itel);
    Itel_dec   = min(Ndec,Itel);
    
    FlagValid(:,Itel) = celestial.scheduling.validate_coo(NewHA(:),NewDec(:),'Lon',InPar.Lon,...
                                                                             'Lat',InPar.Lat,...
                                                                             'AzAltLimit',InPar.AllTelAzAltLimit{Itel_azalt},...
                                                                             'RangeHA',InPar.AllTelRangeHA{Itel_ha},...
                                                                             'RangeDec',InPar.AllTelRangeDec{Itel_dec},...
                                                                             'InterpMethod',InPar.InterpMethod);
end

     


Dist = celestial.coo.sphere_dist_fast(CurRA(:),CurDec(:),NewRA(:).',NewDec(:).');
Dist(~FlagValid) = Inf;

[SD,SI]=sort(Dist.*RAD)

