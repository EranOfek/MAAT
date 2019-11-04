function [Ek,dEdX]=dEdX_calc(Material,Plot)
% Calculate dE/dX as a function of energy
% Package: ultrasat
% Description: Calculate the energy loss of electrons propagating in a
%              material. This is used in order to estimate the Cernekov
%              background.
% Input  : - Material: 'silica' | 'sapphire'. Default is 'silica'.
%          - Plot. Default is true.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - Energy vector [MeV].
%          - dEdX as a function of energy [MeV/(g cm^-2)]
% License: GNU general public license version 3
%     By : Yossi Shvartzvald               Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:  
% Reliable: 2

if nargin<2
    Plot = true;
    if nargin<1
        Material = 'si02_suprasil_2a'; %'silica';   % sapphire
    end
end


%global Silica

e=4.8e-10;
me=9.1e-28;
mp=1.7e-24;
c=3e10;
hbar=1.05e-27;
g=1+10 .^[-3:.01:4];
b=sqrt(1-1 ./g.^2);
Ek=(g-1).*me.*c^2./1.6e-6;


% Silicon
Z=14;
A= 28;
Iav=1.3* 10*Z  *1.6e-12;   % Si
%wp=sqrt(4*pi*rho/A/mp*f0*Z*e^2/me);
dEdXI=2.*pi.*e.^4.*(Z./A)./(mp.*me.*c.^2)./b.^2 /1.6e-6 .* ...
    (log(((g.^2-1).*me.*c.^2./Iav).^2./2 ./(1+g)) - ...
    (2 ./g-1 ./g.^2).*log(2)+1 ./g.^2+(1-1 ./g).^2./8);

dEdXI=max(dEdXI,0);
%dEdXP=2*pi*e^4*(f0*Z/A)/(mp*me*c^2)./b.^2 /1.6e-6 .* ...
%    (log(((g.^2-1)*me*c^2/(wp*hbar)).^2/2 ./(1+g)) - ...
%    (2 ./g-1 ./g.^2)*log(2)+1 ./g.^2+(1-1 ./g).^2/8);
%dEdXP=max(dEdXP,0);
dEdXB=4*Z^2/A*e^6/me^2/mp/c^4/hbar*Ek./b/c*(log(183/Z^(1/3))+1/8);
dEdX=dEdXI+dEdXB; %+dEdXP;

dEdXSi=dEdX;


% Oxygen
Z=8;
A= 16;
Iav=1.3* 10*Z  *1.6e-12;   % O
dEdXI=2*pi*e^4*(Z/A)/(mp*me*c^2)./b.^2 /1.6e-6 .* ...
    (log(((g.^2-1)*me*c^2/Iav).^2/2 ./(1+g)) - ...
    (2 ./g-1 ./g.^2)*log(2)+1 ./g.^2+(1-1 ./g).^2/8);
dEdXI=max(dEdXI,0);
dEdXB=4*Z^2/A*e^6/me^2/mp/c^4/hbar*Ek./b/c*(log(183/Z^(1/3))+1/8);
dEdX=dEdXI+dEdXB; %+dEdXP;

dEdXO=dEdX;

% Aluminium
Z=13;
A= 27;
Iav=1.3* 10*Z  *1.6e-12;   % Al

dEdXI=2*pi*e^4*(Z/A)/(mp*me*c^2)./b.^2 /1.6e-6 .* ...
    (log(((g.^2-1)*me*c^2/Iav).^2/2 ./(1+g)) - ...
    (2 ./g-1 ./g.^2)*log(2)+1 ./g.^2+(1-1 ./g).^2/8);
dEdXI=max(dEdXI,0);
dEdXB=4*Z^2/A*e^6/me^2/mp/c^4/hbar*Ek./b/c*(log(183/Z^(1/3))+1/8);
dEdX=dEdXI+dEdXB; %+dEdXP;

dEdXAl=dEdX;


switch lower(Material)
    case {'sio2','silica','si02_suprasil_2a'}
        %if Silica %SiO2
        dEdX=28/(2*16+28)*dEdXSi+2*16/(2*16+28)*dEdXO;
    case 'sapphire'
        %else %Sapphire Al2O3
        dEdX=27*2/(2*27+3*16)*dEdXAl+3*16/(2*27+3*16)*dEdXO;
        
    otherwise
        error('Unknown Material option');
end


if (Plot)
    figure;
    loglog(Ek,dEdXO);
    hold on;
    loglog(Ek,dEdXSi);
    loglog(Ek,dEdXAl);
    loglog(Ek,dEdX);
    legend('dEdXO','dEdXSi','dEdXAl','dEdX')
    grid
    xlabel('E [MeV]');
    ylabel('dE/dX [MeV/(g cm^{-2})]');


    figure
    loglog(Ek,1./dEdX)
    grid
    axis([1e-2 1e2 1e-2 1])
    xlabel('E [MeV]')
    ylabel('\rho L/E= (dE/dX)^{-1} [(g cm^{-2})/MeV]')
    lw=1.5;
    l_fs=15; 
    sb=12;
    t_fs=15;
    f_w='normal';
    set(gca,'FontName','Times','FontSize',t_fs,'LineWidth',lw,'FontWeight',f_w);
    set(get(gca,'YLabel'),'FontName','Times','FontSize',l_fs,'FontWeight',f_w);
    set(get(gca,'XLabel'),'FontName','Times','FontSize',l_fs,'FontWeight',f_w);
end
