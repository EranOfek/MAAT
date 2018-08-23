function varargout=proper_motion(EpochOut,EpochInRA,EpochInDec,RA,Dec,PM_RA,PM_Dec,Parallax,RadVel)
% Applay proper motion to a catalog
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

Def.Parallax = 1e-4;
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


[Rdot,R] = celestial.coo.pm2space_motion(RA,Dec,PM_RA,PM_Dec,Parallax,RadVel);

if (SameEpoch)
    % assume EpochInRA and EpochInDec are the same
    Rn       = R + bsxfun(@times,Rdot,(EpochOut - EpochInRA));
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
    RnRA        = R + bsxfun(@times,Rdot,(EpochOut - EpochInRA));
    RnDec       = R + bsxfun(@times,Rdot,(EpochOut - EpochInDec));
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


    