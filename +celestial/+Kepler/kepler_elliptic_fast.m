function E=kepler_elliptic_fast(M,Ecc,Tol)
% Solve Kepler equatin (fast version)
% Package: celestial.Kepler
% Description: Fast solver for the Kepler equation for elliptic orbit.
%              This is a simpler version of kepler_elliptic.m
% Input  : - Vector of Mean anomaly [radians].
%          - Vecor of Eccentric anomaly [radians].
%          - Tolerance (e.g., default is 1e-8).
% Output : - Eccentric anomaly.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% see alos: kepler_elliptic.m
% Example: M=rand(1e6,1).*2.*pi; e=rand(1e6,1);
%          tic;E=celestial.Kepler.kepler_elliptic_fast(M,e,1e-8);toc
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3),
    Tol = 1e-8;
end

M = mod(M,2.*pi);

E0 = pi.*ones(size(M)); %M;

%E1 = E0;
E1 = E0 + (M + Ecc.*sin(E0) - E0)./(1 - Ecc.*cos(E0));

%Ns = length(E0);
%IndUnsolved = (1:1:Ns).';

%I = 0;
%Flag = true;
while (any(abs(E1-E0)>Tol)),
    %I = I + 1;
    E0 = E1;
    E1 = E0 + (M + Ecc.*sin(E0) - E0)./(1 - Ecc.*cos(E0));
end
E = E1;

%Diff = E - (M + Ecc.*sin(E));
%max(abs(Diff))

