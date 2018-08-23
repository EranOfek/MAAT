function [X,Y,Z]=moon_elp82(JD_TT,CooType)
% ELP2000-82 ecliptic coordinates of the Moon
% Package: celestial.SolarSys
% Description: Calculate accurate ELP2000-82 ecliptic coordinates of the
%              Moon, referred to the inertial mean ecliptic and equinox of
%              date. This function was previously called moonpos.m.
% Input  : - Vector of julian days in TDT time scale.
%          - Otput coordinates type:
%            'q2000' : equatorial rectangular coordinates referred
%                      to the FK5 equator and equinox
%                      (i.e. mean equator and rotational mean
%                       equinox of J2000).   - (default).
%            'qdate' : equatorial rectangular coordinates referred
%                      to the true equator and equinox of date.
%            'e2000' : ecliptic rectangular coordinates referred to
%                      the inertial mean ecliptic and equinox J2000.0
%            'elp82' : ecliptic rectangular coordinates referred to
%                      the inertial mean ecliptic of date and
%                      departure point \gamma_{2000}'.
%            'elpdt' : ecliptic rectangular coordinates referred to
%                      the inertial mean ecliptic and equinox of date.
% Output : - X [km].
%          - Y [km].
%          - Z [km].
% See Also: mooncool.m
% Reference : ELP2000-82B (Chapront-Touze & Chapront 1982).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jun 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [X,Y,Z]=celestial.SolarSys.moon_elp82(2451545+(1:1:100).');
% Reliable: 2
%--------------------------------------------------------------------------

dms2as = @(D,M,S) D.*3600.0 + M.*60.0 + S;


if (nargin==1),
   CooType = 'q2000';
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end

SizeJD = size(JD_TT);
% rotate JD vector to raw vector.
if (SizeJD(1)==1 && SizeJD(2)==1),
   % do nothing
elseif (SizeJD(1)==1 && SizeJD(2)>1),
   % do nothing
elseif (SizeJD(1)>1 && SizeJD(2)==1),
   JD_TT = JD_TT.';
else
   error('Illegal JD vector size');
end
SizeJD = size(JD_TT);


RAD = 180./pi;
Pc  = 5029.0966;   %["/cy] precession constant in J2000.0

% rotate JD vector to raw vector.

% convert JD_TT to JD_TDB
JD_TDB = JD_TT + celestial.time.tdb_tdt(JD_TT)./86400.0;

% Julian centuary
T = (JD_TDB - 2451545.0)./36525.0;

W1 = dms2as(218,18,59.95571) + 1732559343.73604.*T ...,
                                       - 5.8883.*T.^2 ...,
                                     + 0.006604.*T.^3 ...,
                                     - 0.00003169.*T.^4;
                                  
W1_l = dms2as(218,18,59.95571) + 1732559343.73604.*T;                                  
                                  
W2 = dms2as(83,21,11.67475) + 14643420.2632.*T ...,
                                  - 38.2776.*T.^2 ...,    
                                 - 0.045047.*T.^3 ...,
                               + 0.00021301.*T.^4;
                              
W3 = dms2as(125,2,40.39816) - 6967919.3622.*T ...,
                                  + 6.3622.*T.^2 ...,
                                + 0.007625.*T.^3 ...,
                              - 0.00003586.*T.^4;
                           
% T in original paper                           
M = dms2as(100,27,59.22059) + 129597742.2758.*T ...,
                                    - 0.0202.*T.^2 ...,
                                  + 0.000009.*T.^3 ...,
                                + 0.00000015.*T.^4;
                             
M_l = dms2as(100,27,59.22059) + 129597742.2758.*T;                             
                             
Om = dms2as(102,56,14.42753) + 1161.2283.*T ...,
                                + 0.5327.*T.^2 ...,
                              - 0.000138.*T.^3;

D = dms2as(297,51,0.73512) + 1602961601.4603.*T ...,
                                    - 5.8681.*T.^2 ...,
                                  + 0.006595.*T.^3 ...,
                                - 0.00003184.*T.^4;
                             
D_l = dms2as(297,51,0.73512) + 1602961601.4603.*T;
                                                                
Lt = dms2as(357,31,44.79306) + 129596581.0474.*T ...,
                                     - 0.5529.*T.^2 ...,
                                     + 0.000147.*T.^3;
                                  
Lt_l = dms2as(357,31,44.79306) + 129596581.0474.*T;                             
                                
L = dms2as(134,57,48.28096) + 1717915923.4728.*T ...,
                                    + 32.3893.*T.^2 ...,
                                   + 0.051651.*T.^3 ...,
                                   - 0.00024470.*T.^4;
                                
L_l = dms2as(134,57,48.28096) + 1717915923.4728.*T;
                              
F = dms2as(93,16,19.55755) + 1739527263.0983.*T ...,
                                   - 12.2505.*T.^2 ...,
                                  - 0.001021.*T.^3 ...,
                                + 0.00000417.*T.^4;
                             
F_l = dms2as(93,16,19.55755) + 1739527263.0983.*T;                            

Zeta = W1_l + Pc.*T;

% planetary longitude (VSOP82)
Me   = dms2as(252,15,3.25986)  + 538101628.68898.*T;
V    = dms2as(181,58,47.28305) + 210664136.43355.*T;
Ma   = dms2as(355,25,59.78866) +  68905077.59284.*T;
J    = dms2as(34,21,5.34212)   +  10925660.42861.*T;
S    = dms2as(50,4,38.89694)   +   4399609.65932.*T;
U    = dms2as(314,3,18.01841)  +   1542481.19393.*T;
N    = dms2as(304,20,55.19575) +    786550.32074.*T;



% convert arcsec to radians
W1   = W1./(3600.*RAD);
W2   = W2./(3600.*RAD);
W3   = W3./(3600.*RAD);
M    = M./(3600.*RAD);
M_l  = M_l./(3600.*RAD);
Om   = Om./(3600.*RAD);
D    = D./(3600.*RAD);
D_l  = D_l./(3600.*RAD);
Lt   = Lt./(3600.*RAD);
Lt_l = Lt_l./(3600.*RAD);
L    = L./(3600.*RAD);
L_l  = L_l./(3600.*RAD);
F    = F./(3600.*RAD);
F_l  = F_l./(3600.*RAD);
Zeta = Zeta./(3600.*RAD);
Me   = Me./(3600.*RAD);
V    = V./(3600.*RAD);
Ma   = Ma./(3600.*RAD);
J    = J./(3600.*RAD);
S    = S./(3600.*RAD);
U    = U./(3600.*RAD);
N    = N./(3600.*RAD);


%D  = W1 - M + pi;
%L  = W1 - W2;
%F  = W1 - W3;

% converting to [0, 2.*pi] range
con2range = @(X) (X./(2.*pi) - floor(X./(2.*pi))).*2.*pi;
W1   = con2range(W1);
W2   = con2range(W2);
W3   = con2range(W3);
M    = con2range(M);
M_l  = con2range(M_l);
Om   = con2range(Om);
D    = con2range(D);
D_l  = con2range(D_l);
Lt   = con2range(Lt);
Lt_l = con2range(Lt_l);
L    = con2range(L);
L_l  = con2range(L_l);
F    = con2range(F);
F_l  = con2range(F_l);
Zeta = con2range(Zeta);
Me   = con2range(Me);
V    = con2range(V);
J    = con2range(J);
S    = con2range(S);
U    = con2range(U);
N    = con2range(N);

OnesJD = ones(SizeJD);
T2Vec  = T.^2;


Nu  = 1732559343.18;   %"/cy
Nt  =  129597742.34;   %"/cy
Ms  = Nt./Nu;          % m
Det = -0.12879;        % "
DE  =  0.01789;        % "
DG  = -0.08066;        % "
Dnt = -0.0642;         % "/cy
Dnu =  0.55604;        % "/cy
Alp =  0.0026;   %???



% load files
load ELP1.mat
load ELP2.mat
load ELP3.mat

%ELP1-3

DA1 = -Ms.*(ELP1(:,6) + 2.*Alp.*ELP1(:,10)./(3.*Ms)).*(Dnu./Nu) + ...,
           (ELP1(:,6) + 2.*Alp.*ELP1(:,10)./(3.*Ms)).*(Dnt./Nu) + ...,
       (ELP1(:,7).*DG + ELP1(:,8).*DE + ELP1(:,9).*Det)./206264.81;
DA2 = -Ms.*(ELP2(:,6) + 2.*Alp.*ELP2(:,10)./(3.*Ms)).*(Dnu./Nu) + ...,
           (ELP2(:,6) + 2.*Alp.*ELP2(:,10)./(3.*Ms)).*(Dnt./Nu) + ...,
       (ELP2(:,7).*DG + ELP2(:,8).*DE + ELP2(:,9).*Det)./206264.81;
DA3 = -Ms.*(ELP3(:,6) + 2.*Alp.*ELP3(:,10)./(3.*Ms) + 2.*ELP3(:,5)./(3.*Ms)).*Dnu/Nu + ...,
           (ELP3(:,6) + 2.*Alp.*ELP3(:,10)./(3.*Ms)).*Dnt./Nu + ...,
       (ELP3(:,7).*DG + ELP3(:,8).*DE + ELP3(:,9).*Det)./206264.81;
    
Sl = sum(((ELP1(:,5)+DA1)*OnesJD).*sin(ELP1(:,1)*D + ELP1(:,2)*Lt + ELP1(:,3)*L + ELP1(:,4)*F));
Sb = sum(((ELP2(:,5)+DA2)*OnesJD).*sin(ELP2(:,1)*D + ELP2(:,2)*Lt + ELP2(:,3)*L + ELP2(:,4)*F));
Sr = sum(((ELP3(:,5)+DA3)*OnesJD).*cos(ELP3(:,1)*D + ELP3(:,2)*Lt + ELP3(:,3)*L + ELP3(:,4)*F));

load ELP4.mat
load ELP5.mat
load ELP6.mat

%ELP4-6
Sl = Sl + sum((ELP4(:,7)*OnesJD).*sin(ELP4(:,1)*Zeta + ELP4(:,2)*D_l + ELP4(:,3)*Lt_l + ELP4(:,4)*L_l + ELP4(:,5)*F_l + (ELP4(:,6)*OnesJD)./RAD));
Sb = Sb + sum((ELP5(:,7)*OnesJD).*sin(ELP5(:,1)*Zeta + ELP5(:,2)*D_l + ELP5(:,3)*Lt_l + ELP5(:,4)*L_l + ELP5(:,5)*F_l + (ELP5(:,6)*OnesJD)./RAD));
Sr = Sr + sum((ELP6(:,7)*OnesJD).*sin(ELP6(:,1)*Zeta + ELP6(:,2)*D_l + ELP6(:,3)*Lt_l + ELP6(:,4)*L_l + ELP6(:,5)*F_l + (ELP6(:,6)*OnesJD)./RAD));

load ELP7.mat
load ELP8.mat
load ELP9.mat

%ELP7-9
Sl = Sl + T.*sum((ELP7(:,7)*OnesJD).*sin(ELP7(:,1)*Zeta + ELP7(:,2)*D_l + ELP7(:,3)*Lt_l + ELP7(:,4)*L_l + ELP7(:,5)*F_l + (ELP7(:,6)*OnesJD)./RAD));
Sb = Sb + T.*sum((ELP8(:,7)*OnesJD).*sin(ELP8(:,1)*Zeta + ELP8(:,2)*D_l + ELP8(:,3)*Lt_l + ELP8(:,4)*L_l + ELP8(:,5)*F_l + (ELP8(:,6)*OnesJD)./RAD));
Sr = Sr + T.*sum((ELP9(:,7)*OnesJD).*sin(ELP9(:,1)*Zeta + ELP9(:,2)*D_l + ELP9(:,3)*Lt_l + ELP9(:,4)*L_l + ELP9(:,5)*F_l + (ELP9(:,6)*OnesJD)./RAD));

load ELP10.mat
load ELP11.mat
load ELP12.mat

%ELP10-12
Sl = Sl + sum((ELP10(:,13)*OnesJD).*sin(ELP10(:,1)*Me + ELP10(:,2)*V +ELP10(:,3)*M_l + ELP10(:,4)*Ma + ELP10(:,5)*J + ELP10(:,6)*S + ELP10(:,7)*U + ELP10(:,8)*N + ELP10(:,9)*D_l + ELP10(:,10)*L_l + ELP10(:,11)*F_l + (ELP10(:,12)*OnesJD)./RAD));
Sb = Sb + sum((ELP11(:,13)*OnesJD).*sin(ELP11(:,1)*Me + ELP11(:,2)*V +ELP11(:,3)*M_l + ELP11(:,4)*Ma + ELP11(:,5)*J + ELP11(:,6)*S + ELP11(:,7)*U + ELP11(:,8)*N + ELP11(:,9)*D_l + ELP11(:,10)*L_l + ELP11(:,11)*F_l + (ELP11(:,12)*OnesJD)./RAD));
Sr = Sr + sum((ELP12(:,13)*OnesJD).*sin(ELP12(:,1)*Me + ELP12(:,2)*V +ELP12(:,3)*M_l + ELP12(:,4)*Ma + ELP12(:,5)*J + ELP12(:,6)*S + ELP12(:,7)*U + ELP12(:,8)*N + ELP12(:,9)*D_l + ELP12(:,10)*L_l + ELP12(:,11)*F_l + (ELP12(:,12)*OnesJD)./RAD));

load ELP13.mat
load ELP14.mat
load ELP15.mat

%ELP13-15
Sl = Sl + T.*sum((ELP13(:,13)*OnesJD).*sin(ELP13(:,1)*Me + ELP13(:,2)*V +ELP13(:,3)*M_l + ELP13(:,4)*Ma + ELP13(:,5)*J + ELP13(:,6)*S + ELP13(:,7)*U + ELP13(:,8)*N + ELP13(:,9)*D_l + ELP13(:,10)*L_l + ELP13(:,11)*F_l + (ELP13(:,12)*OnesJD)./RAD));
Sb = Sb + T.*sum((ELP14(:,13)*OnesJD).*sin(ELP14(:,1)*Me + ELP14(:,2)*V +ELP14(:,3)*M_l + ELP14(:,4)*Ma + ELP14(:,5)*J + ELP14(:,6)*S + ELP14(:,7)*U + ELP14(:,8)*N + ELP14(:,9)*D_l + ELP14(:,10)*L_l + ELP14(:,11)*F_l + (ELP14(:,12)*OnesJD)./RAD));
Sr = Sr + T.*sum((ELP15(:,13)*OnesJD).*sin(ELP15(:,1)*Me + ELP15(:,2)*V +ELP15(:,3)*M_l + ELP15(:,4)*Ma + ELP15(:,5)*J + ELP15(:,6)*S + ELP15(:,7)*U + ELP15(:,8)*N + ELP15(:,9)*D_l + ELP15(:,10)*L_l + ELP15(:,11)*F_l + (ELP15(:,12)*OnesJD)./RAD));

load ELP16.mat
load ELP17.mat
load ELP18.mat

%ELP16-18
Sl = Sl + sum((ELP16(:,13)*OnesJD).*sin(ELP16(:,1)*Me + ELP16(:,2)*V +ELP16(:,3)*M_l + ELP16(:,4)*Ma + ELP16(:,5)*J + ELP16(:,6)*S + ELP16(:,7)*U + ELP16(:,8)*D_l + ELP16(:,9)*Lt_l + ELP16(:,10)*L_l + ELP16(:,11)*F_l + (ELP16(:,12)*OnesJD)./RAD));
Sb = Sb + sum((ELP17(:,13)*OnesJD).*sin(ELP17(:,1)*Me + ELP17(:,2)*V +ELP17(:,3)*M_l + ELP17(:,4)*Ma + ELP17(:,5)*J + ELP17(:,6)*S + ELP17(:,7)*U + ELP17(:,8)*D_l + ELP17(:,9)*Lt_l + ELP17(:,10)*L_l + ELP17(:,11)*F_l + (ELP17(:,12)*OnesJD)./RAD));
Sr = Sr + sum((ELP18(:,13)*OnesJD).*sin(ELP18(:,1)*Me + ELP18(:,2)*V +ELP18(:,3)*M_l + ELP18(:,4)*Ma + ELP18(:,5)*J + ELP18(:,6)*S + ELP18(:,7)*U + ELP18(:,8)*D_l + ELP18(:,9)*Lt_l + ELP18(:,10)*L_l + ELP18(:,11)*F_l + (ELP18(:,12)*OnesJD)./RAD));

load ELP19.mat
load ELP20.mat
load ELP21.mat

%ELP19-21
Sl = Sl + T.*sum((ELP19(:,13)*OnesJD).*sin(ELP19(:,1)*Me + ELP19(:,2)*V +ELP19(:,3)*M_l + ELP19(:,4)*Ma + ELP19(:,5)*J + ELP19(:,6)*S + ELP19(:,7)*U + ELP19(:,8)*D_l + ELP19(:,9)*Lt_l + ELP19(:,10)*L_l + ELP19(:,11)*F_l + (ELP19(:,12)*OnesJD)./RAD));
Sb = Sb + T.*sum((ELP20(:,13)*OnesJD).*sin(ELP20(:,1)*Me + ELP20(:,2)*V +ELP20(:,3)*M_l + ELP20(:,4)*Ma + ELP20(:,5)*J + ELP20(:,6)*S + ELP20(:,7)*U + ELP20(:,8)*D_l + ELP20(:,9)*Lt_l + ELP20(:,10)*L_l + ELP20(:,11)*F_l + (ELP20(:,12)*OnesJD)./RAD));
Sr = Sr + T.*sum((ELP21(:,13)*OnesJD).*sin(ELP21(:,1)*Me + ELP21(:,2)*V +ELP21(:,3)*M_l + ELP21(:,4)*Ma + ELP21(:,5)*J + ELP21(:,6)*S + ELP21(:,7)*U + ELP21(:,8)*D_l + ELP21(:,9)*Lt_l + ELP21(:,10)*L_l + ELP21(:,11)*F_l + (ELP21(:,12)*OnesJD)./RAD));

load ELP22.mat
load ELP23.mat
load ELP24.mat

%ELP22-24
Sl = Sl + sum((ELP22(:,7)*OnesJD).*sin(ELP22(:,2)*D_l + ELP22(:,3)*Lt_l + ELP22(:,4)*L_l + ELP22(:,5)*F_l + (ELP22(:,6)*OnesJD)./RAD));
Sb = Sb + sum((ELP23(:,7)*OnesJD).*sin(ELP23(:,2)*D_l + ELP23(:,3)*Lt_l + ELP23(:,4)*L_l + ELP23(:,5)*F_l + (ELP23(:,6)*OnesJD)./RAD));
Sr = Sr + sum((ELP24(:,7)*OnesJD).*sin(ELP24(:,2)*D_l + ELP24(:,3)*Lt_l + ELP24(:,4)*L_l + ELP24(:,5)*F_l + (ELP24(:,6)*OnesJD)./RAD));

load ELP25.mat
load ELP26.mat
load ELP27.mat

%ELP25-27
Sl = Sl + T.*sum((ELP25(:,7)*OnesJD).*sin(ELP25(:,2)*D_l + ELP25(:,3)*Lt_l + ELP25(:,4)*L_l + ELP25(:,5)*F_l + (ELP25(:,6)*OnesJD)./RAD));
Sb = Sb + T.*sum((ELP26(:,7)*OnesJD).*sin(ELP26(:,2)*D_l + ELP26(:,3)*Lt_l + ELP26(:,4)*L_l + ELP26(:,5)*F_l + (ELP26(:,6)*OnesJD)./RAD));
Sr = Sr + T.*sum((ELP27(:,7)*OnesJD).*sin(ELP27(:,2)*D_l + ELP27(:,3)*Lt_l + ELP27(:,4)*L_l + ELP27(:,5)*F_l + (ELP27(:,6)*OnesJD)./RAD));

load ELP28.mat
load ELP29.mat
load ELP30.mat

%ELP22-24
Sl = Sl + sum((ELP28(:,7)*OnesJD).*sin(ELP28(:,2)*D_l + ELP28(:,3)*Lt_l + ELP28(:,4)*L_l + ELP28(:,5)*F_l + (ELP28(:,6)*OnesJD)./RAD));
Sb = Sb + sum((ELP29(:,7)*OnesJD).*sin(ELP29(:,2)*D_l + ELP29(:,3)*Lt_l + ELP29(:,4)*L_l + ELP29(:,5)*F_l + (ELP29(:,6)*OnesJD)./RAD));
Sr = Sr + sum((ELP30(:,7)*OnesJD).*sin(ELP30(:,2)*D_l + ELP30(:,3)*Lt_l + ELP30(:,4)*L_l + ELP30(:,5)*F_l + (ELP30(:,6)*OnesJD)./RAD));

load ELP31.mat
load ELP32.mat
load ELP33.mat

%ELP22-24
Sl = Sl + sum((ELP31(:,7)*OnesJD).*sin(ELP31(:,2)*D_l + ELP31(:,3)*Lt_l + ELP31(:,4)*L_l + ELP31(:,5)*F_l + (ELP31(:,6)*OnesJD)./RAD));
Sb = Sb + sum((ELP32(:,7)*OnesJD).*sin(ELP32(:,2)*D_l + ELP32(:,3)*Lt_l + ELP32(:,4)*L_l + ELP32(:,5)*F_l + (ELP32(:,6)*OnesJD)./RAD));
Sr = Sr + sum((ELP33(:,7)*OnesJD).*sin(ELP33(:,2)*D_l + ELP33(:,3)*Lt_l + ELP33(:,4)*L_l + ELP33(:,5)*F_l + (ELP33(:,6)*OnesJD)./RAD));

load ELP34.mat
load ELP35.mat
load ELP36.mat

%ELP34-36
Sl = Sl + T2Vec.*sum((ELP34(:,7)*OnesJD).*sin(ELP34(:,2)*D_l + ELP34(:,3)*Lt_l + ELP34(:,4)*L_l + ELP34(:,5)*F_l + (ELP34(:,6)*OnesJD)./RAD));
Sb = Sb + T2Vec.*sum((ELP35(:,7)*OnesJD).*sin(ELP35(:,2)*D_l + ELP35(:,3)*Lt_l + ELP35(:,4)*L_l + ELP35(:,5)*F_l + (ELP35(:,6)*OnesJD)./RAD));
Sr = Sr + T2Vec.*sum((ELP36(:,7)*OnesJD).*sin(ELP36(:,2)*D_l + ELP36(:,3)*Lt_l + ELP36(:,4)*L_l + ELP36(:,5)*F_l + (ELP36(:,6)*OnesJD)./RAD));

% convert arcsec to radians
switch CooType
 case 'elpdt'
    %MoonLon = W1 + (Sl + Pc.*T)./(3600.*RAD);
    MoonLon = W1 + Sl./(3600.*RAD);
 otherwise
    MoonLon = W1 + Sl./(3600.*RAD);
end
MoonLat = Sb./(3600.*RAD);
MoonRad = Sr;

% coordinates are in the inertial mean ecliptic of date and departure point \gamma_{2000}^{'}


% convert to range [0, 2*pi]
MoonLon = con2range(MoonLon);
MoonLat = con2range(MoonLat);


switch CooType
 case {'q2000','J2000'}
    %  'q2000' : equatorial rectangular coordinates referred
    %            to the FK5 equator and equinox
    %            (i.e. mean equator and rotational mean
    %            equinox of J2000).
    [Xe,Ye,Ze] = elp2e2000(T,MoonLon,MoonLat,MoonRad);
    [X,Y,Z]    = e2000q2000(Xe,Ye,Ze);
 case {'qdate','date'}
    %  'qdate' : equatorial rectangular coordinates referred
    %            to the true equator and equinox of date
    [Xe,Ye,Ze] = elp2e2000(T,MoonLon,MoonLat,MoonRad);
    [X,Y,Z]    = e2000q2000(Xe,Ye,Ze);
    for I=1:1:SizeJD,
       RotM = rotm_coo('p',JD_TT(I));
       XYZ  = RotM*[X(I);Y(I);Z(I)];
       X(I) = XYZ(1);
       Y(I) = XYZ(2);
       Z(I) = XYZ(3);
    end
       
 case {'e2000','eJ2000'}
    %  'e2000' : ecliptic rectangular coordinates referred to
    %            the inertial mean ecliptic and equinox J2000.0
    [X,Y,Z] = elp2e2000(T,MoonLon,MoonLat,MoonRad);

 case {'elp82'}  
    %  'elp82' : ecliptic rectangular coordinates referred to
    %            the inertial mean ecliptic of date and
    %            departure point \gamma_{2000}'.
    X = MoonRad.*cos(MoonLon).*cos(MoonLat);
    Y = MoonRad.*sin(MoonLon).*cos(MoonLat);
    Z = MoonRad.*sin(MoonLat);
    
 case {'elpdt','edate'}   
    %  'elpdt' : ecliptic rectangular coordinates referred to
    %            the inertial mean ecliptic and equinox of date.
    [X,Y,Z] = elp2date(T,MoonLon,MoonLat,MoonRad);
   
 otherwise
    error('Illegal CooType');
end






%--------------------------------------------
% coordinates conversion section
%--------------------------------------------

%-----------------------------------------------------
function [X,Y,Z]=elp2date(T,MoonLon,MoonLat,MoonRad)
% convert to inertial mean ecliptic of date and equinox of date.
RAD = 180./pi;

Pa = 5029.0966.*T + 1.1120.*T.^2 + 0.000077.*T.^3 - 0.00002353.*T.^4;
Pa = Pa./(3600.*RAD);

MoonLon = MoonLon + Pa;

MoonLon = 2.*pi.*(MoonLon./(2.*pi) - floor(MoonLon./(2.*pi)));

X = MoonRad.*cos(MoonLon).*cos(MoonLat);
Y = MoonRad.*sin(MoonLon).*cos(MoonLat);
Z = MoonRad.*sin(MoonLat);



%-----------------------------------------------------
function [X,Y,Z]=elp2e2000(T,MoonLon,MoonLat,MoonRad)
% convert to inertial mean ecliptic of date and equinox of date.

Xe = MoonRad.*cos(MoonLon).*cos(MoonLat);
Ye = MoonRad.*sin(MoonLon).*cos(MoonLat);
Ze = MoonRad.*sin(MoonLat);

% rotate for eact T
RotM = zeros(3,3);
X = zeros(size(T));
Y = zeros(size(T));
Z = zeros(size(T));
for I=1:1:length(T),
   P =  0.10180391e-4.*T(I) + ...,
        0.47020439e-6.*T(I).^2 - ...,
        0.5417367e-9.*T(I).^3 - ...,
        0.2507948e-11.*T(I).^4 + ...,
        0.463486e-14.*T(I).^5;
   Q = -0.113469002e-3.*T(I) + ...,
        0.12372674e-6.*T(I).^2 + ...,
        0.12654170e-8.*T(I).^3 - ...,
        0.1371808e-11.*T(I).^4 - ...,
        0.320334e-14.*T(I).^5;
  
   RotM(1,1) = 1-2.*P.^2;
   RotM(1,2) = 2.*P.*Q;
   RotM(1,3) = 2.*P.*sqrt(1-P.^2-Q.^2);
   RotM(2,1) = 2.*P.*Q;
   RotM(2,2) = 1-2.*Q.^2;
   RotM(2,3) = -2.*Q.*sqrt(1-P.^2-Q.^2);
   RotM(3,1) = -2.*P.*sqrt(1-P.^2-Q.^2);
   RotM(3,2) = 2.*Q.*sqrt(1-P.^2-Q.^2);
   RotM(3,3) = 1-2.*P.^2-2.*Q.^2;

   PosVec = RotM*[Xe(I);Ye(I);Ze(I)];
   X(I) = PosVec(1);
   Y(I) = PosVec(2);
   Z(I) = PosVec(3);
end

%-----------------------------------------------------
function [X,Y,Z]=e2000q2000(Xe,Ye,Ze)
% convert from inertial mean ecliptic of date and equinox of date
% to FK5 equatorial.

RotM = zeros(3,3);
RotM(1,1) =  1.000000000000;
RotM(1,2) =  0.000000437913;
RotM(1,3) = -0.000000189859;
RotM(2,1) = -0.000000477299;
RotM(2,2) =  0.917482137607;
RotM(2,3) = -0.397776981701;
RotM(3,1) =  0.000000000000;
RotM(3,2) =  0.397776981701;
RotM(3,3) =  0.917482137607;

PosVec = RotM*[Xe;Ye;Ze];
X = PosVec(1,:);
Y = PosVec(2,:);
Z = PosVec(3,:);




%function ArcSec=dms2as(D,M,S)
%-----------------------------------------------------------------
% In arguments  : - Deg.
%          - Min.
%          - Sec.
% Out    : - arcsec.
%ArcSec = D.*3600.0 + M.*60.0 + S;



