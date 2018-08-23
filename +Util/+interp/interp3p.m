function [Fp,ExtramP,ExtramFp,Extram2D]=interp3p(Y,P)
% Stirling interpolation
% Package: Util.interp
% Description: Given three, equally spaced, data points (Y), and a vector
%              of normalized positions in units in which the three data
%              points are at [-1,0,1], return, based on the Stirling
%              interpolation formula, the values of the interpolation
%              function at the vector of positions, and information
%              regarding the extramum of the interpolation function.
% Input  : - Vector of three values specified at positions [-1,0,1].
%          - Vector of positions at which to calculate the values
%            of the interpolation function.
% Output : - Values of the interpolation function.
%          - Position of extramum.
%          - Value at extramum.
%          - 2nd derivative at extramum.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Fp,ExtramP,ExtramFp,Extram2D]=Util.interp.interp3p([10,1,8],[-1:0.01:1]');
% Reliable: 
%--------------------------------------------------------------------------

Dph = Y(3)-Y(2);  % delta_{+1/2}
Dmh = Y(2)-Y(1);  % delta_{-1/2}
D20 = Dph - Dmh;  % delta^{2}_{0}

Par = [0.5.*D20, 0.5.*(Dph+Dmh), Y(2)];
if (~isempty(P))
    Fp  = polyval(Par,P);
else
    Fp = [];
end

ExtramP  = -Par(2)./(2.*Par(1));
ExtramFp = Par(1).*ExtramP.^2 + Par(2).*ExtramP + Par(3); %polyval(Par,ExtramP);
Extram2D = 2.*Par(1);






