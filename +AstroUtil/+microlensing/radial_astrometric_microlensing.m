function [Shift,Res]=radial_astrometric_microlensing(Mass,D_l,D_s,D_ls,Beta0,PM,T,T0)
% The astrometric deflection of the primary microleneing image.
% Package: AstroUtil.microlensing
% Description: Calculate the astrometric deflection of the primary
%              microlensing image as a function of time.
% Input  : - Deflector mass [solar mass].
%          - D_l [pc].
%          - D_s [pc].
%          - D_ls [pc].
%          - Minimum impact parameter at T0: Beta0 [arcsec].
%          - Proper motion [arcsec/unit time].
%          - Vector of times [unit time].
%          - T0 [unit time].
% Output : - Shifts in arcsec as a function of time T.
%          - A structure array with the following fields:
%            'ER' - The Einstein radius [arcsec].
%            'BetaT' - Vector of deflector source distance [arcsec].
%            'T'  - Vector of times.;


% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [S,Res]=AstroUtil.microlensing.radial_astrometric_microlensing
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (nargin==0)
    Mass  = 1./50000; 1;
    D_l   = 1./206; 10;
    D_s   = 1000;
    D_ls  = D_s - D_l;
    Beta0 = 2;   % [arcsec]
    PM    = 1;    % [arcsec/yr]
    T0    = 0;
    T     = (-10:1:10)';
    
end

U2    = Beta0.^2 + (PM.*(T-T0)).^2;
BetaT = sqrt(U2);

%Nt = numel(BetaT);
%for It=1:1:Nt

[ER,T1,T2,Mu1,Mu2,TD1,TD2,Tcm]=AstroUtil.microlensing.pointsource_lens(Mass,D_l,D_s,D_ls,BetaT./(RAD.*3600));
Shift = (T1.*RAD.*3600 - BetaT);
Res.ER    = ER.*RAD.*3600;
Res.BetaT = BetaT;
Res.T     = T;
