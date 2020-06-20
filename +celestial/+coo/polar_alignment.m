function [Res]=polar_alignment(Meas_H,Meas_D,Measured_Dt_dt,Phi,Plot)
% Calculate the RA/Dec drift due to equatorial polar alignemnt error.
% Package: celestial
% Description: 
% Input  : - Vector of HA of targets [rad].
%          - Vector Dec of targets [rad].
%          - Vector of measured Declination drifts ["/s].
%          - True altitude of NCP (Geodetic latitude of observer) [rad].
% Output : - Structure of results.
%            .BestAlt - Alt shift in [deg] needed to fix polar alignment.
%                   Positive upward.
%            .BestAz - Az shift in [deg] needed to fix polar alignment
%                   Positive Eastward.
%            .BestPsi - Best fit hour angle of telescope pole relative to
%                   NCP [deg], measured westward from the superior merdian.
%            .BestBeta - Best fit angular distance [deg] between telescope
%                   pole and NCP.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res]=celestial.coo.polar_alignment
% Reliable: 
%--------------------------------------------------------------------------

if nargout<5
    Plot = false;
end

RAD = 180./pi;


% D - true Dec
% Dt - telescope Dec.
% R - true RA
% Rt - telescope RA
% H - true HA
% Ht - telescope HA
% Phi - true altitude of NCP (latitude)
% Phit - telescope polar altitude
% Z - zenith distance of star
% Beta - distance between NCP and telescope pole
% Psi - HA of telescope pole

if nargin==0
    % simulation mode
    VecHA = (-90:10:90);
    VecDec = (-20:10:80);
    [MatH,MatD] = meshgrid(VecHA./RAD,VecDec./RAD);
    Phi = 32./RAD;
    Psi = -100./RAD;
    Beta = 0.55./RAD;

    [dHt_dt,dDt_dt,Ht,Dt]=celestial.coo.polar_alignment_drift(MatH,MatD,Phi,Psi,Beta);
    %Drift = sqrt(dRt_dt.^2 + dDt_dt.^2);
    %surface(VecHA,VecDec,dDt_dt)

    % simulate measurments with noise
    RelNoise = 0.01;
    Meas_H = [-60:30:60]'; %[-75:25:75].'; %[-75:25:75]';
    Meas_D = 0.*ones(size(Meas_H));
    Measured_Dt_dt = interp2(VecHA,VecDec,dDt_dt,Meas_H,Meas_D);
    Measured_Dt_dt = Measured_Dt_dt +randn(size(Measured_Dt_dt)).*Measured_Dt_dt.*RelNoise;

end

    
D_Beta  = [0.1:0.01:2]'./RAD;
D_Psi   = [-180:2:180]'./RAD;
Nbeta   = numel(D_Beta);
Npsi    = numel(D_Psi);
RMS     = zeros(Nbeta,Npsi);
for Ibeta=1:1:Nbeta
    for Ipsi=1:1:Npsi
        [~,Pred_dDt_dt]=celestial.coo.polar_alignment_drift(Meas_H./RAD,Meas_D./RAD,Phi,D_Psi(Ipsi),D_Beta(Ibeta));
        
        RMS(Ibeta,Ipsi) = std(Measured_Dt_dt - Pred_dDt_dt);
    end
end

%Plot = true
if Plot
    surface(D_Psi.*RAD,D_Beta.*RAD,RMS);
    colorbar
end
[Min,MinI]=Util.stat.minnd(RMS);
Res.BestBeta = D_Beta(MinI(1));
Res.BestPsi  = D_Psi(MinI(2));
Res.Phit     = asin( cos(Res.BestBeta).*sin(Phi) + sin(Res.BestBeta).*cos(Phi).*cos(Res.BestPsi) );
Res.BestDAlt = -(Res.Phit - Phi);
Res.BestDAz  = asin(sin(Res.BestBeta).*sin(Res.BestPsi)./cos(Res.Phit));

Res.BestBeta = Res.BestBeta.*RAD;
Res.BestPsi  = Res.BestPsi.*RAD;
Res.BestDAlt = Res.BestDAlt.*RAD;
Res.BestDAz  = Res.BestDAz.*RAD;

% 
% %% old
% 
% Phit = asin(sin(Phi)*cos(Beta) + cos(Phi).*sin(Beta).*cos(Psi));
% Alt = celestial.coo.ha2alt(H,D,Phi);
% Z   = pi./2 - Alt;
% 
% dH_dt = 360.*3600./86164.091; % ["/s]   sidereal rate
% dH_dt = dH_dt./(3600.*RAD);
% 
% Dt = asin(cos(Beta).*sin(D) + sin(Beta).*cos(D).*cos(Psi - H));
% 
% dDt_dt = dH_dt .* (cos(D).*sin(Beta).*sin(Psi-H))./cos(Dt);
% 
% Ht = acos(cos(Z)./(cos(Dt).*cos(Phit)) - tan(Dt).*tan(Phit));
% 
% % note that there is an error in the last term of Equation 4 in Markworth
% dHt_dt = dH_dt .* cos(D).*cos(Phi).*sin(H)./(cos(Dt).*cos(Phit).*sin(Ht)) - sec(Dt).^2.*dDt_dt.*cos(Z).*sin(Dt)./(sin(Ht).*cos(Phit)) + ...
%          dDt_dt.*tan(Phit)./(cos(Dt).^2 .* sin(Ht));
% 
%      
% dDt_dt = dDt_dt.*3600.*RAD;
% dRt_dt = (dHt_dt - dH_dt).*3600.*RAD;
