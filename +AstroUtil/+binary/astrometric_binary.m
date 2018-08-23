function [Orbit,Par]=astrometric_binary(JD,El)
% Predict astrometric binary position
% Package: AstroUtil.binary
% Description: Given orbital elements of an astrometric binary, predicts
%              its sky position as a function of time.
% Input  : - JD or dates (see julday.m for options) in which to calculate
%            the binary position.
%          - Binary orbital elements. This is either a vector of
%            [Period (days), a (arcsec), e, T (JD),...
%                            Omega (rad),omega (rad) ,i (rad)],
%            or a structure containing the following fields:
%            .Period, .a, .e, .T, .Omega, .omega, .i.
%            Default is [365 1 0.5 1 1 1 1].
% Output : - Structue containing the orbit position as a function of time.
%            The following fields are available:
%            .JD - JD
%            .X, .Y, .Z - positions.
%          - Structure containing the following additional parameters:
%            .Nu   - true anomaly
%            .R    - radius vector
%            .E    - Eccentric anomaly
%            .dNudt- d\nu/dt
% Tested : Maylab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Orbit,Par]=AstroUtil.binary.astrometric_binary; plot(Orbit.X,Orbit.Y,'.');
% Reliable: 2
%--------------------------------------------------------------------------

Def.JD = 2451545 + (0:1:365).';
Def.El = [365 1 0.5 1 1 1 1];

if (nargin==0)
    JD = Def.JD;
    El = Def.El;
elseif (nargin==1)
    El = Def.El;
elseif (nargin==2)
    % do nothing
else
    error('Illegal number of input arguments');
end
    
if (size(JD,2)==1)
    % do nothing
    Orbit.JD = JD;
else
   Orbit.JD = celestial.time.julday(JD);
end

if (isnumeric(El))
    ElS.Period = El(1);
    ElS.a      = El(2);
    ElS.e      = El(3);
    ElS.T      = El(4);
    ElS.Omega  = El(5);
    ElS.omega  = El(6);
    ElS.i      = El(7);
else
    ElS = El;
end


ElS.n   = 2.*pi./ElS.Period;  % [rad/day]
M       = ElS.n.*(Orbit.JD - ElS.T);
M       = celestial.coo.angle_in2pi(M);
ElS.q   = ElS.a.*(1 - ElS.e);
[Par.Nu,Par.R,Par.E]  = celestial.Kepler.kepler_elliptic(M, ElS.q, ElS.e, NaN);
Par.dNudt = ElS.n.*sqrt(1-ElS.e.^2)./((1-ElS.e.*cos(Par.E)).^2);
%[Nu,R,dNudt]
[Orbit.X,Orbit.Y,Orbit.Z] = celestial.Kepler.trueanom2pos(Par.R,Par.Nu,ElS.Omega,ElS.omega,ElS.i);


