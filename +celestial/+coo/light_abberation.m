function P2=light_abberation(P1,Et)
%--------------------------------------------------------------------------
% light_abberation function                                          ephem
% Description: Given an object observer-centric direction, corrected
%              for light deflection in the natural frame (P1),
%              calculate the proper direction of the object (P2) in
%              the observer-centric inertial frame that is moving
%              with instantaneous velocity of the observer (V)
%              relative to the natural frame.
% Input  : - (P1) object observer-centric direction, corrected
%            for light deflection in the natural frame. [au]
%            (3-element column vector).
%            If more than one column is given, then the proper direction
%            is calculated for each vector. In that case Q and E
%            should contain the same number of columns.
%          - (Et) Observer velocity in the barycentric frame. [au/day]
% Output : - The object proper direction (P2) in the
%            observer-centric inertial frame that is moving
%            with instantaneous velocity of the observer (V)
%            relative to the natural frame.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% Reliable:
%--------------------------------------------------------------------------
C = 173.1446327205364; % speed of light [au/day]

V    = Et./C;  % convert to speed of light units
AbsV = ones(3,1)*sqrt(sum((V.*V),1));
Beta = 1./sqrt(1-AbsV.^2);
     
% scalar products
P1Vsp = ones(3,1)*sum(P1.*V,1);
Temp = 1 + (P1Vsp)./(1+1./Beta);

P2 = (P1./Beta + Temp.*V)./(1+P1Vsp);

P2 = P2./[ones(3,1)*sqrt(sum(P2.*P2))];


