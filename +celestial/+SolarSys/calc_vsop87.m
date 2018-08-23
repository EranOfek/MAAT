function [Coo,Vel]=calc_vsop87(JD_TT, ObjectName, VSOP_Type, OutCooType)
% Planetary coordinates based on the VSOP87 theory
% Package: celestial.SolarSys
% Description: Calculate planetary coordinates using the VSOP87 theory.
% Input  : - Column vector of JD in the TT time scale.
%          - Object name:
%            'Mercury' : Mercury
%            'Venus'   : Venus
%            'Earth'   : Earth
%            'EMB'     : Earth-Moon barycenter
%            'Mars'    : Mars
%            'Jupiter' : Jupiter
%            'Saturn'  : Saturn
%            'Uranus'  : Uranus
%            'Neptune' : Neptune
%            'Sun'     : Sun
%          - Type of VSOP87 series to use, VSOP_Type option:
%            'a' : Heliocentric; ec. rectangular; equinox and ecliptic of J2000.0
%            'b' : Heliocentric; ec. spherical; equinox and ecliptic of J2000.0
%            'c' : Heliocentric; ec. rectangular; equinox and ecliptic of date
%            'd' : Heliocentric; ec. spherical; equinox and ecliptic of date
%            'e' : Barycentric; ec. rectangular; equinox and ecliptic of J2000.0
%          - Type of output coordinates:
%            'd' : default coordinates as given by the series type (default).
%            'E' : rectangular, Equatorial FK5 J2000.0 coordinate system,
%                  possible in the case of VSOP_Type='A'|'B'|'E'.
% Output : - Coordinates matrix, column per time.
%            In each column [X;Y;Z] or [L;B;R]. X,Y,Z,R in au. L,B in radinas.
%          - Velocity matrix, column per time.
%            In each column [Xt;Yt;Zt] or [Lt;Bt;Rt]. Xt,Yt,Zt,Rt in au/day.
%            Lt,Bt in radinas/day.
% Referebce: Bretagnon P., Francou G. 1988 A&A 202, 309
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Coo,Vel]=celestial.SolarSys.calc_vsop87(2451545+(0:1:100)', 'Uranus', 'e', 'E');
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin==3),
   OutCooType = 'd';
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end



switch VSOP_Type
 case 'a'
    CooType = 'R';  % rectangular
    System  = 'J';  % J2000.0
 case 'c'
    CooType = 'R';  % rectangular
    System  = 'D';  % date
 case 'e'
    CooType = 'R';  % rectangular
    System  = 'J';  % J2000.0
 case 'b'
    CooType = 'S';  % spherical
    System  = 'J';  % J2000.0
 case 'd'
    CooType = 'S';  % spherical
    System  = 'D';  % date
 otherwise
    error('Unknwon VSOP87 series type');
end
 
 
 
FileBaseName = 'vsop87';
switch lower(ObjectName)
 case 'mercury'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.mer.mat'];
 case 'venus'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.ven.mat'];
 case 'earth'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.ear.mat'];    
 case 'emr'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.emb.mat'];
 case 'mars'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.mar.mat'];
 case 'jupiter'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.jup.mat'];
 case 'saturn'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.sat.mat'];
 case 'uranus'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.ura.mat'];
 case 'neptune'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.nep.mat'];
 case 'sun'
    VSOP87_FileName = [FileBaseName,VSOP_Type,'.sun.mat'];
 otherwise
    error('Unknown object name');
end
    


TJY = 365250;  % 1000 Julian years
LOD = 86400;   % length of day [sec].
T = (JD_TT - 2451545.0)/TJY;
Dt   = 1./(LOD.*TJY);
Dts  = 1./(LOD);
Tdt1 = T - Dt;
Tdt2 = T + Dt;


N = length(T);

% rotation matrix for transforming rectangular coordinates
% defined in dynamical ecliptic J2000.0 frame (i.e., VSOP87A, VSOP87E)
% to equatorial frame FK5 J2000.0
Rot_Dyn2FK5 = [+1.000000000000  +0.000000440360  -0.000000190919;
               -0.000000479966  +0.917482137087  -0.397776982902;
               +0.000000000000  +0.397776982902  +0.917482137087];



Lambda = zeros(N,12);
Lambda(:,1)  = 4.40260884240 + 26087.9031415742 .* T;  %  Mercury
Lambda(:,2)  = 3.17614669689 + 10213.2855462110 .* T;  %    Venus
Lambda(:,3)  = 1.75347045953 +  6283.0758499914 .* T;  %    Earth
Lambda(:,4)  = 6.20347611291 +  3340.6124266998 .* T;  %    Mars
Lambda(:,5)  = 0.59954649739 +   529.6909650946 .* T;  %    Jupiter
Lambda(:,6)  = 0.87401675650 +   213.2990954380 .* T;  %    Saturn
Lambda(:,7)  = 5.48129387159 +    74.7815985673 .* T;  %    Uranus
Lambda(:,8)  = 5.31188628676 +    38.1330356378 .* T;  %    Neptune
Lambda(:,9)  = 5.19846674103 + 77713.7714681205 .* T;  %    Moon D
Lambda(:,10) = 1.62790523337 + 84334.6615813083 .* T;  %   Moon F
Lambda(:,11) = 2.35555589827 + 83286.9142695536 .* T;  %   Moon l
Lambda(:,12) = 3.81034454697 + 83997.0911355954 .* T;  %   Moon Lm


S_Col = 15;
K_Col = 16;
A_Col = 17;
B_Col = 18;
C_Col = 19;
P_Col = [3:14];
% read modified VSOP87 file
%VSOP87 = load(VSOP87_FileName);
load(VSOP87_FileName);
eval(sprintf('VSOP87 = %s%s;',FileBaseName,VSOP_Type));

Pow = mod(VSOP87(:,1),10);
Var = floor(mod(VSOP87(:,1),100)./10);
Obj = floor(mod(VSOP87(:,1),1000)./100);
Typ = floor(VSOP87(:,1)./1000);

%Coo1    = zeros(3,N);
Coo     = zeros(3,N);
Vel     = zeros(3,N);
Coo_dt1 = zeros(3,N);
Coo_dt2 = zeros(3,N);

for CooInd=1:1:3,
   J = find(Var==CooInd);
   for I=1:1:N,
      % calc coordinates per time
      Coo(CooInd,I) = sum((T(I).^(Pow(J))).*VSOP87(J,A_Col).*cos(VSOP87(J,B_Col) + VSOP87(J,C_Col).*T(I)));

      %Phi = sum(VSOP87(J,P_Col)*(Lambda(I,:).'),2);
      %Coo1(CooInd,I) = sum((T(I).^(Pow(J))).*(VSOP87(J,S_Col).*sin(Phi) + VSOP87(J,K_Col).*cos(Phi)));
      % calc coordinates for time+dt (for velocity)
      Coo_dt1(CooInd,I) = sum((Tdt1(I).^(Pow(J))).*VSOP87(J,A_Col).*cos(VSOP87(J,B_Col) + VSOP87(J,C_Col).*Tdt1(I)));
      Coo_dt2(CooInd,I) = sum((Tdt2(I).^(Pow(J))).*VSOP87(J,A_Col).*cos(VSOP87(J,B_Col) + VSOP87(J,C_Col).*Tdt2(I)));
   end
end

Vel = 0.5.*(Coo_dt2 - Coo_dt1)./Dts;
if (CooType=='S'),
   Coo(1:2,:) = (Coo(1:2,:)./(2.*pi) - floor(Coo(1:2,:)./(2.*pi))).*2.*pi;
   Ip = find(Coo(2,:)>pi);
   Coo(2,Ip) = Coo(2,Ip) - 2.*pi;
end


switch OutCooType
 case 'd'
    % --- do nothing ---
 case 'E'
    % --- convert to Equatorial FK5 J2000.0 ---
    % rotation matrix for transforming rectangular coordinates
    % defined in dynamical ecliptic J2000.0 frame (i.e., VSOP87A, VSOP87E)
    % to equatorial frame FK5 J2000.0
    Rot_Dyn2FK5 = [+1.000000000000  +0.000000440360  -0.000000190919;
                   -0.000000479966  +0.917482137087  -0.397776982902;
                   +0.000000000000  +0.397776982902  +0.917482137087];
    % rotate coordinates
    switch CooType
     case 'R'
        Coo = Rot_Dyn2FK5*Coo;
        Vel = Rot_Dyn2FK5*Vel;
     case 'S'
        X  = cos(Coo(2,:)).*cos(Coo(1,:));
        Y  = cos(Coo(2,:)).*sin(Coo(1,:));
        Z  =                sin(Coo(1,:));
        Xt = cos(Coo(2,:)).*cos(Coo(1,:));
        Yt = cos(Coo(2,:)).*sin(Coo(1,:));
        Zt =                sin(Coo(1,:));        
        Coo = Rot_Dyn2FK5*[X;Y;Z];
        Vel = Rot_Dyn2FK5*[Xt;Yt;Zt];
     otherwise
        error ('Unknown CooType');
     end
 otherwise('Unknown OutCooType option');
end
            
