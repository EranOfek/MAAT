function jup_satcurve(Month,Year)
% Plot monthly curves of the position of the Galilean satellites
% Package: celestial.SolarSys
% Description: Plot monthly curves of the relative position of the
%              Galilean satellites of Jupiter.
% Input  : - Month
%          - Year
% Output : null
% plot   : - Jovian Satellites curve for one month.
% Reference : Meeus, J. 1991 in: Astronomical Algorithms.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Web example: http://astroclub.tau.ac.il/ephem/JovSat/
% Example: celestial.SolarSys.jup_satcurve(1,2015)
% Reliable: 1
%--------------------------------------------------------------------------
TimeStep  = 0.01;
JD_Start  = celestial.time.julday([0, Month, Year, 0]);
JD        = (JD_Start:TimeStep:JD_Start+32)';
Time      = (0:TimeStep:32)';
JupRad    = 0.5;

Color1 = 'r';
Color2 = 'g';
Color3 = 'b';
Color4 = 'm';

RAD = 180./pi;

d   = JD - 2451545.0;
V   = (172.74  + 0.00111588.*d)./RAD;
M   = (357.529 + 0.9856003.*d)./RAD;
N   = (20.020  + 0.0830853.*d + 0.329.*sin(V))./RAD;
J   = (66.115  + 0.9025179.*d - 0.329.*sin(V))./RAD;
A   = (1.915.*sin(M) + 0.020.*sin(2.*M))./RAD;
B   = (5.555.*sin(N) + 0.168.*sin(2.*N))./RAD;
K   = J + A - B;
R   = 1.00014 - 0.01671.*cos(M) - 0.00014.*cos(2.*M);
r   = 5.20872 - 0.25208.*cos(N) - 0.00611.*cos(2.*N);
Del = sqrt(r.^2 + R.^2 - 2.*r.*R.*cos(K));
Psi = asin(R.*sin(K)./Del);


u1  = (163.8067 + 203.4058643.*(d - Del./173))./RAD + Psi - B;
u2  = (358.4108 + 101.2916334.*(d - Del./173))./RAD + Psi - B;
u3  = (  5.7129 +  50.2345179.*(d - Del./173))./RAD + Psi - B;
u4  = (224.8151 +  21.4879801.*(d - Del./173))./RAD + Psi - B;

G   = (331.18 + 50.310482.*(d - Del./173))./RAD;
H   = ( 87.40 + 21.569231.*(d - Del./173))./RAD;

r1  =  5.9073 - 0.0244.*cos(2.*(u1-u2));
r2  =  9.3991 - 0.0882.*cos(2.*(u2-u3));
r3  = 14.9924 - 0.0216.*cos(G);
r4  = 26.3699 - 0.1935.*cos(H);

u1  = u1 + 0.473.*sin(2.*(u1-u2))./RAD;
u2  = u2 + 1.065.*sin(2.*(u2-u3))./RAD;
u3  = u3 + 0.165.*sin(G)./RAD;
u4  = u4 + 0.841.*sin(H)./RAD;

Lam   = (34.35 + 0.083091.*d + 0.329.*sin(V))./RAD + B;
Ds    = 3.12.*sin(Lam + 42.8./RAD);
De    = Ds - 2.22.*sin(Psi).*cos(Lam+22./RAD) - 1.30.*(r-Del).*sin(Lam-100.5./RAD)./Del;
Ds    = Ds./RAD;
De    = De./RAD;

X1  =  r1.*sin(u1);
Y1  = -r1.*cos(u1).*sin(De);
X2  =  r2.*sin(u2);
Y2  = -r2.*cos(u2).*sin(De);
X3  =  r3.*sin(u3);
Y3  = -r3.*cos(u3).*sin(De);
X4  =  r4.*sin(u4);
Y4  = -r4.*cos(u4).*sin(De);

% motion direction
DX1 = diff(X1);
DX1 = [DX1;DX1(length(DX1))];
DX2 = diff(X2);
DX2 = [DX2;DX2(length(DX2))];
DX3 = diff(X3);
DX3 = [DX3;DX3(length(DX3))];
DX4 = diff(X4);
DX4 = [DX4;DX4(length(DX4))];

Nt    = length(Time);
Ht    = floor(Nt./2);
OffSet = 35;
Htime  = 16.0;
Time1 = rot90(Time(1:Ht),2);
Time2 = rot90(Time(Ht+1:Nt),2);
Boundary = 70;

axis('off');
hold on;

X1L = X1(1:Ht);
X2L = X2(1:Ht);
X3L = X3(1:Ht);
X4L = X4(1:Ht);

X1R = X1(Ht+1:Nt);
X2R = X2(Ht+1:Nt);
X3R = X3(Ht+1:Nt);
X4R = X4(Ht+1:Nt);

DX1L = DX1(1:Ht);
DX2L = DX2(1:Ht);
DX3L = DX3(1:Ht);
DX4L = DX4(1:Ht);

DX1R = DX1(Ht+1:Nt);
DX2R = DX2(Ht+1:Nt);
DX3R = DX3(Ht+1:Nt);
DX4R = DX4(Ht+1:Nt);


set(gcf,'Color','White');

axis([-Boundary Boundary 0 Htime+3]);

% plot title
TitleStr = ['Satellites of Jupiter, ',num2str(Month),'/',num2str(Year)];
h_t = text(-Boundary./2,Htime+1.5,TitleStr);
set(h_t,'FontSize',12);
h_t = text(Boundary-20,Htime+1.5,'By: Eran O. Ofek');
set(h_t,'FontSize',5);

% plot left chart
I1 = find(~(DX1L<0 & abs(X1L)<JupRad));
h_p = plot(-X1L(I1)-OffSet,Time1(I1),'r.','Color',Color1);
set(h_p,'MarkerSize',1);

I2 = find(~(DX2L<0 & abs(X2L)<JupRad));
h_p = plot(-X2L(I2)-OffSet,Time1(I2),'g.','Color',Color2);
set(h_p,'MarkerSize',1);

I3 = find(~(DX3L<0 & abs(X3L)<JupRad));
h_p = plot(-X3L(I3)-OffSet,Time1(I3),'b.','Color',Color3);
set(h_p,'MarkerSize',1);

I4 = find(~(DX4L<0 & abs(X4L)<JupRad));
h_p = plot(-X4L(I4)-OffSet,Time1(I4),'y.','Color',Color4);
set(h_p,'MarkerSize',1);

% plot right chart
I1 = find(~(DX1R<0 & abs(X1R)<JupRad));
h_p = plot(-X1R(I1)+OffSet,Time2(I1)-Htime,'r.','Color',Color1);
set(h_p,'MarkerSize',1);

I2 = find(~(DX2R<0 & abs(X2R)<JupRad));
h_p = plot(-X2R(I2)+OffSet,Time2(I2)-Htime,'g.','Color',Color2);
set(h_p,'MarkerSize',1);

I3 = find(~(DX3R<0 & abs(X3R)<JupRad));
h_p = plot(-X3R(I3)+OffSet,Time2(I3)-Htime,'b.','Color',Color3);
set(h_p,'MarkerSize',1);

I4 = find(~(DX4R<0 & abs(X4R)<JupRad));
h_p = plot(-X4R(I4)+OffSet,Time2(I4)-Htime,'y.','Color',Color4);
set(h_p,'MarkerSize',1);

% plot jupiter
plot([-OffSet-JupRad;-OffSet-JupRad],[0;Htime],'k');
plot([-OffSet+JupRad;-OffSet+JupRad],[0;Htime],'k');
plot([OffSet-JupRad;OffSet-JupRad],[0;Htime],'k');
plot([OffSet+JupRad;OffSet+JupRad],[0;Htime],'k');

% plot boundary lines
plot([-3;-3],[0;Htime+0.5],'k');
plot([Boundary-3;Boundary-3],[0;Htime+0.5],'k');
h_t = text(-Boundary+5,Htime+0.3,'West');
set(h_t,'FontSize',7);
h_t = text(5,Htime+0.3,'West');
set(h_t,'FontSize',7);
h_t = text(-12,Htime+0.3,'East');
set(h_t,'FontSize',7);
h_t = text(Boundary-12,Htime+0.3,'East');
set(h_t,'FontSize',7);



% plot days lines and labels
Day = 0;
for Dy=Htime:-1:(Htime-floor(Htime))
   plot([-Boundary+3;-3],[Dy;Dy],':k');
   plot([3;Boundary-3],[Dy;Dy],':k');
   h_t = text(-Boundary,Dy,num2str(Day));
   set(h_t,'FontSize',6);
   h_t = text(0,Dy,num2str(Day+floor(Htime)));
   set(h_t,'FontSize',6);
   Day = Day + 1;
end   

hold off;
