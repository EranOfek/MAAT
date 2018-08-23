function [N,EdgesS]=lognlogs(S,EdgesS);
%---------------------------------------------------------------------------
% lognlogs function                                               AstroStat
% Description: Calculate Log N - Log S plot for set of flux measurments.
% Input  : - List of flux measurments.
%          - List of Edges in which to calculate N(>S).
% Output : - Cumulative number of fluexs above S [e.g., N(>S)].
%          - Vector of fluxes.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                   January 2008
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: D=rand(100000,3); D=sqrt(sum(D.^2,2)); I=find(D<1);
%          D=D(I); L=D.^-2;
%          [N,S] = lognlogs(L,logspace(1,4,30)');
%          loglog(S,N);
%          Par=fitpow(S,N,sqrt(N));
%---------------------------------------------------------------------------

NE = length(EdgesS);
N  = zeros(NE,1);
for I=1:1:NE,
   Ks = find(S>EdgesS(I));
   N(I) = length(Ks);
end
