function [DAz,DAlt,Best_d,Best_Ep]=polar_alignment(h,Dec,Beta,DeltaT)
% SHORT DESCRIPTION HERE
% Package: celestial
% Description: 
% Input  : - HA [deg].
%          - Dec [deg]
%          - Beta [arcsec]
%          - DeltaT [seconds]
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DAz,DAlt,d,Ep]=celestial.coo.polar_alignment
%          [DAz,DAlt,d,Ep]=celestial.coo.polar_alignment
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

if nargin==0
%     % simulation mode
%     h = [-90:30:90].';
%     Dec = zeros(size(h));
%     d   = 0.3
%     Ep  = 30;
%     gamma = acot((cosd(Dec).*cotd(d) - sind(Dec).*cosd(Ep+h))./sind(Ep+h));  % rad
%     DeltaT = 60;  % sec
%     Beta = sin(gamma).*DeltaT./240;
  

    % plot mode
    d = 1;
    Ep = 0;
    
    h = [-180:0.5:180].';
    D = [-62:0.5:85].';
    [Math,MatD]=meshgrid(h,Dec);
    gamma = acot((cosd(MatD).*cotd(d) - sind(MatD).*cosd(Ep+Math))./sind(Ep+Math));  % rad
    Alt = celestial.coo.ha2alt(Math./RAD,MatD./RAD,30./RAD);
    Az  = celestial.coo.ha2az(Math./RAD,MatD./RAD,30./RAD);
    
    axesm ('eqaazim', 'Frame', 'on', 'Grid', 'on');
    Flag = Alt>0;
    pcolorm(Alt(Flag).*RAD,Az(Flag).*RAD,gamma(Flag))
    
end




    

% Beta in arcseconds
% DeltaT in seconds
gamma = asin(Beta./(15.*DeltaT));  % radians

% solve:
% cot(gamma) = (cos(delta) cot(d) - sin(delta) cos(epsilon+h))/sin(epslion+h)


Vec_d  = (0.01:0.01:3);
Vec_Ep = (10:5:360);
Nd     = numel(Vec_d);
Nep    = numel(Vec_Ep);
Resid  = zeros(Nd,Nep);
for Id=1:1:Nd
    ResidI = cot(gamma) - (cosd(Dec).*cotd(Vec_d(Id)) - sind(Dec).*cosd(Vec_Ep+h))./sind(Vec_Ep+h);
    ResidI(ResidI==Inf | ResidI==-Inf) = NaN;
    Resid(Id,:) = sqrt(nansum(ResidI.^2,1));
end
[Min,MinI]=Util.stat.minnd(Resid);
Best_d  = Vec_d(MinI(1));
Best_Ep = Vec_Ep(MinI(2));

DAz  = Best_d.*sind(Best_Ep);
DAlt = Best_d.*cosd(Best_Ep);

contour(Vec_Ep,Vec_d,log10(Resid));
H= xlabel('$\epsilon$ [deg]');
H.FontSize = 18;
H.Interpreter = 'latex';
H= ylabel('$d$ [deg]');
H.FontSize = 18;
H.Interpreter = 'latex';

