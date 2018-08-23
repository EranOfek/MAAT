function [Nu,R]=kepler_lowecc(M,e,Q)
% A low eccentricity serise solution for the Kepler equation
% Package: celestial.Kepler
% Description: Solve the Kepler Equation for low eccentricity using
%              a series approximation. This is typically better than 1"
%              for e<0.1.
% Input  : - Mean anomaly [rad].
%          - Eccentricity (e).
%          - Periastron distance (q).
% Output : - True anomaly [rad].
% Reference: Meeus (1991)
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Nu]=celestial.Kepler.kepler_lowecc(1,0.01)
%          [Nu,R]=celestial.Kepler.kepler_lowecc(1,0.01,1); 
% Reliable: 2
%--------------------------------------------------------------------------

% Nu = M + 2.*e.*sin(M) + 5./4.*e.^2.*sin(2.*M) + ...
%     e.^3./12 .*(13.*sin(3.*M)-3.*sin(M)) +...
%     e.^4./96.*(103.*sin(4.*M) - 44.*sin(2.*M));

Nu = M + (2.*e - 0.25.*e.^3 + 5./96.*e.^5).*sin(M) + ...
         (1.25.*e.^2 - 11./24.*e.^4).*sin(2.*M) + ...
         (13./12.*e.^3 - 43./64.*e.^5).*sin(3.*M) + ...
         (103./96.*e.^4).*sin(4.*M) + ...
         (1097./960.*e.^5).*sin(5.*M);

if (nargout>1),     
    R = Q.*(1+e)./(1+e.*cos(Nu));
end
