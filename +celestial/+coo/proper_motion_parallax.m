function varargout=proper_motion_parallax(EpochOut,EpochInRA,EpochInDec,RA,Dec,PM_RA,PM_Dec,Parallax,RadVel)
% Applay proper motion and parallaxto a catalog
% Package: celestial.coo
% Description: Applay proper motion to a catalog
% Input  : - Final epoch (days). Either a scalar or a vector.
%          - Catalog initial epoch (days) of RA.
%            Either a scalar or a vector.
%          - Catalog initial epoch (days) of Dec. If empty then will
%            use the EpochIn of RA.
%            Either a scalar or a vector.
%          - RA at initial epoch [radians].
%          - Dec at initial epoch [radians].
%          - Proper motion in RA [mas/yr].
%          - Proper motion in Dec [mas/yr].
%          - Parallax [mas], default is 1e-4;
%          - Radial velocity [km/s], default is 0.
% Output * Either 2 or 3 output arguments.
%          If two output arguments, these are the final RA and Dec
%          in radians.
%          [RA,Dec] = proper_motion(...);
%          If three output arguments, these are the X,Y,Z cosine directions
%          in units of AU.
%          [X,Y,Z] = proper_motion(...);
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jan 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------
arsc2rad=pi/180/3600;
Def.Parallax = 0;
Def.RadVel   = 0;
if (nargin==7)
   Parallax = Def.Parallax;
   RadVel   = Def.RadVel;
elseif (nargin==8)
   RadVel   = Def.RadVel;
elseif (nargin==9)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(EpochInDec))
    EpochInDec = EpochInRA;
    SameEpoch  = true;
else
    SameEpoch  = false;
end

% Rdot [rad/day] , R [cosined unity vector]
[Rdot,R] = pm2space_motion_rad(RA,Dec,PM_RA,PM_Dec,Parallax,RadVel);

if (SameEpoch)
    % assume EpochInRA and EpochInDec are the same
    %Calculate the bayrcentric position for the Out epoch
    [Coo,Vel] = celestial.SolarSys.calc_vsop87(EpochOut, 'Earth', 'e', 'E');
    %Calculate the bayrcentric position for the in epoch (catalog epoch)
    [CooRef,Vel] = celestial.SolarSys.calc_vsop87(EpochInRA(1), 'Earth', 'e', 'E');
    
    % pi Eb [mas au] - the correction of the barycentric position at the
    % output epoch (piEb) and catalog epoch 
    
    piEb=     Parallax .* Coo' ; %masc * au
    piEbRef=     Parallax .* CooRef' ; %masc * au
    % pi Eb [rad au]
    piEb =    piEb* arsc2rad/1000;
    piEbRef = piEbRef* arsc2rad/1000;
    Rn       = R + bsxfun(@times,Rdot,(EpochOut - EpochInRA))-piEb +piEbRef ;
    if (nargout==2)
       [varargout{1},varargout{2}] = celestial.coo.cosined2coo(Rn(:,1),Rn(:,2),Rn(:,3));
    elseif (nargout==3)
       varargout{1} = Rn(:,1);
       varargout{2} = Rn(:,2);
       varargout{3} = Rn(:,3);
    else
        error('Illegal number of output arguments');
    end
else
    % assume EpochInRA and EpochInDec are different
    %Calculate the bayrcentric position for the Out epoch
    [Coo,Vel] = celestial.SolarSys.calc_vsop87(EpochOut, 'Earth', 'e', 'E');
    
    %Calculate the bayrcentric position for the in epoch (catalog epoch RA and Dec)
    [CooRefRA,Vel] = celestial.SolarSys.calc_vsop87(EpochInRA(1), 'Earth', 'e', 'E');
    [CooRefDec,Vel] = celestial.SolarSys.calc_vsop87(EpochInDec(1), 'Earth', 'e', 'E');
    %{
    if (numel(unique(EpochInRA))~=1)
        %rotate the eliptict coordinates to be equotorial
        for i= 1:numel(EpochInRA)
            RotM=celestial.coo.rotm_coo('E',EpochInRA(i));
            CooRA(:,i)= RotM*CooRA(:,i);
            CooDec(:,i)= RotM*CooDec(:,i);
        end
    end
    %} 
    
    % pi Eb [mas au] - the correction of the barycentric position at the
    % output epoch (piEb) and catalog epoch RA and Dec separetly
    piEb=    Parallax .* Coo'; 
    piEbRefRA = Parallax .* CooRefRA'; 
    piEbRefDec = Parallax .* CooRefDec'; 
    
    % pi Eb - >[rad au]
    piEb=    piEb* arsc2rad/1000;
    piEbRefRA = piEbRefRA* arsc2rad/1000;
    piEbRefDec = piEbRefDec* arsc2rad/1000;
    
    RnRA        = R + bsxfun(@times,Rdot,(EpochOut - EpochInRA))    -piEb + piEbRefRA;
    RnDec       = R + bsxfun(@times,Rdot,(EpochOut - EpochInDec))   -piEb + piEbRefDec;
    
    
    [OutRA,Tmp] = celestial.coo.cosined2coo(RnRA(:,1),RnRA(:,2),RnRA(:,3));    
    [~,OutDec]  = celestial.coo.cosined2coo(RnDec(:,1),RnDec(:,2),RnDec(:,3));
    if (nargout==2)
        varargout{1} = OutRA;
        varargout{2} = OutDec;
     elseif (nargout==3)
       [varargout{1}, varargout{2}, varargout{3}] = celestial.coo.coo2cosined(OutRA,OutDec);
    else
        error('Illegal number of output arguments');   
    end   
end



function [SpaceMotion,SpaceVec]=pm2space_motion_rad(RA,Dec,PM_RA,PM_Dec,Parallax,RadVel)
%--------------------------------------------------------------------------
% pm2space_motion function                                           ephem
% Description: Convert proper motion, radial velocity and parralax to
%              space motion vector in the equatorial system.
% Input  : - [RA] in radians.
%          - [Dec] in radians.
%          - PM in RA [mas/yr]
%          - PM in Dec [mas/yr]
%          - Parallax [mas]
%          - Radial velocity [km/s]
% Output : - Space motion vector [X Y Z] in au/day, in the equatorial
%            system.
%          - Space position vector in au, in the equatorial system. 
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Seidelmann 1992 p.121
% Reliable: 2
%--------------------------------------------------------------------------
RAD  = 180./pi;
S    = 2.*pi./(360.*3600.*36525);          % convert from "arcsec/cy to radians/day
AU   = constant.au('SI')./1000;      % [km]
K    = 86400./AU;                          % convert from km/s to au/day

% convert PM_RA from mas/yr to second of arcsec/cy
PM_RA  = PM_RA .* 100./1000 ./cos(Dec);


% conver PM_Dec from mas/yr to arcsec/cy
PM_Dec = PM_Dec.* 100./1000;

% make sure parallax is not zero
Parallax(abs(Parallax)<1e-5) = 1e-5;   % mas

% Convert Parallax from mas to radians
Parallax = Parallax ./(1000.*3600.*RAD);

N = max([length(RA);length(Dec);length(PM_RA);length(PM_Dec);length(Parallax);length(RadVel)]);
RA      = RA.*ones(N,1);
Dec     = Dec.*ones(N,1);
PM_RA   = PM_RA.*ones(N,1);
PM_Dec  = PM_Dec.*ones(N,1);
Parallax= Parallax.*ones(N,1);
RadVel  = RadVel.*ones(N,1);

SpaceVec = celestial.coo.cosined([RA,Dec]);

SpaceMotion = zeros(N,3);
for I=1:1:N

   Mat = [-cos(Dec(I)).*sin(RA(I)), -sin(Dec(I)).*cos(RA(I)), cos(Dec(I)).*cos(RA(I)); ...
           cos(Dec(I)).*cos(RA(I)), -sin(Dec(I)).*sin(RA(I)), cos(Dec(I)).*sin(RA(I)); ...
	   0                      ,  cos(Dec(I))            , sin(Dec(I))];

   Vec = [S.*PM_RA(I); S.*PM_Dec(I); Parallax(I).*K.*RadVel(I)];

   SpaceMotion(I,:) = (Mat*Vec).';
end


    