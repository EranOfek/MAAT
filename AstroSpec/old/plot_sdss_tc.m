function plot_sdss_tc;
%----------------------------------------------------------------------
% plot_sdss_tc function                                      AstroSpec
% Description: Plot SDSS filters transmission curves
% Input  : null
% Output : null
% Plot   : SDSS ugriz filters transmission curves
% Tested : Matlab 5.3
%     By : Eran O. Ofek         Feb. 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: plot_sdss_tc;
%----------------------------------------------------------------------

U = get_filter('SDSS','u');
G = get_filter('SDSS','g');
R = get_filter('SDSS','r');
I = get_filter('SDSS','i');
Z = get_filter('SDSS','z');


plot(U.nT{1}(:,1),U.nT{1}(:,2),'m');
hold on;
plot(G.nT{1}(:,1),G.nT{1}(:,2),'b');
plot(R.nT{1}(:,1),R.nT{1}(:,2),'r');
plot(I.nT{1}(:,1),I.nT{1}(:,2),'y');
plot(Z.nT{1}(:,1),Z.nT{1}(:,2),'k');
