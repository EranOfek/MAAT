function [FlagIn,H]=plot_polygonselect(X,Y,varargin)
%--------------------------------------------------------------------------
% plot_polygonselect function                                     plotting
% Description: Plot and let the use
% Input  : - 
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    May 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


H = plot(X,Y,varargin{:});


B=1;
Ip = 0;
while (B==1)
    Ip = Ip+1;
    [XB(Ip),YB(Ip),B] = ginput(1);
    if (Ip>1),
        hold on;
        plot([XB(Ip-1),XB(Ip)],[YB(Ip-1),YB(Ip)],'r--');
    end
end


FlagIn = inpolygon(X,Y,XB,YB);
