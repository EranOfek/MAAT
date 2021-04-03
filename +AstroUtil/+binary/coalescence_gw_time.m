function [t_coal,da_dt,de_dt]=coalescence_gw_time(M1,M2,a,e)
% Calculate the coalescence time for a binary stars due to GW emission
% Package: +AstroUtil.binary
% Input  : - M1 [Solar mass]
%          - M2 [Solar mass]
%          - a [cm]
%          - e. Default is 0. (Correct for low e only).
% Output : - Coalescensc time [s]
%          - da/dt [cm/s]
%          - de/dt [s^-1]
% Example: [t_coal,da_dt,de_dt]=AstroUtil.binary.coalescence_gw_time(M1,M2,a,e)

if nargin<4
    e = 0;
end

M1 = M1.*constant.SunM;
M2 = M2.*constant.SunM;


Fe = (1-e.^2).^(-7./2) .* (1 + e.^2 .* 73./24 + e.^4 .* 37./96);
da_dt = -64./5 .* constant.G.^3./(constant.c.^5) .* M1.*M2.*(M1 + M2).*Fe./(a.^3);
de_dt = -304./15 .* constant.G.^3./(constant.c.^5) .* M1.*M2.*(M1+M2)./(a.^4)  .*e.*(1-e.^2).^(-5./2) .* (1 + e.^2 .*121./304);
t_coal = 5./256 .*constant.c.^5./(constant.G.^3).*a.^4./(M1.*M2.*(M1+M2).*Fe);


if (1==0)
    M2=0.6;
    M1=logspace(log10(0.6),log10(30),100)';
    p=logspace(log10(0.1),log10(100),300).*60;

    R=celestial.Kepler.kepler3law(2e33.*(M1+M2),'p',p);
    [t_coal,da_dt,de_dt]=AstroUtil.binary.coalescence_gw_time(M1,M2,R.a,0);
    rH=(M1./(3.*(M1+M2))).^(1/3);
    q=M2./M1;
    rRL=0.49.*q.^(2./3) ./(0.6.*q.^(2./3) + log(1+q.^(1./3)));

    FlagRL=R.a.*rRL>6400e5;     
    FlagRH=R.a.*rH>6400e5;     

    [c,h]=contour(p,M1,log10(t_coal./(365.25.*86400)));
    clabel(c,h);
    set(gca','YS','log','XS','log')
    hold on;
    contour(p,M1,FlagRL,'k-')
    contour(p,M1,FlagRH,'k--')
    H=xlabel('Period [s]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    H=ylabel('M1 [solar mass]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    
end
