function solve_alignment6dof_problem(Info)
%


% DOF:
% Phi   - HA of the mount-axis as measured Eastward from the NCP
% Theta - HA of the telescope-axis as measured Eastward from the axis
% D     - NCP-Axis angular distance
% R     - Telescope-Axis angular distance
% L0    - encoder Dec zero
% H0    - encoder HA zero

% additional parameters
% H     - HA of pointing relative to axis
% Ht    - Ht = H - H0
% L     - pointing-axis angular distance
% Lt    - Lt = L - L0

% Input:
% HA_Tel, Dec_Tel - from astrometry
% HA_Axis, Dec_Axis - from mount

RAD = 180./pi;

load Info_good.mat
Info = Info_good;

N = numel(Info);
for I=1:1:N
    HA_Axis(I)  = Info(I).HA;
    Dec_Axis(I) = Info(I).Dec;
    
    RADec = celestial.coo.coco(Info(I).Asmtry_coo,'j2000.0','J2021.0');
    
    HA_Tel(I)   = celestial.time.lst(Info(I).JD,34.81./RAD).*360 - RADec(1).*RAD;
    Dec_Tel(I)  = RADec(2).*RAD;
end


Flag           = Dec_Axis>90;
Dec_Axis(Flag) = 90-(Dec_Axis(Flag)-90);
HA_Axis(Flag)  = HA_Axis(Flag) + 180;
HA_Axis        = mod(HA_Axis,360);

L = 90 - Dec_Axis;
H = HA_Axis;


Lat_NCP = 90;
Lon_NCP = 0;


if 1==0
    Vec_Phi   = [(-180:30:179).'];
    Vec_Theta = (0:30:359).';
    Vec_D     = (0.0:0.5:2.5).';
    Vec_R     = (0.0:0.5:3.5).';
    Vec_L0    = (-10:1:10).';
    Vec_H0    = (-180:30:179).';



    Vec_Phi   = [(0:2:10).'];
    Vec_Theta = (0:5:359).';
    Vec_D     = (1.5:0.1:2.5).';
    Vec_R     = (0:0.1:3).';
    Vec_L0    = (-2:0.1:4).';
    Vec_H0    = (-10:2:40).';

    Vec_Phi   = (0:2:10).';
    Vec_Theta = (0:2:10).';
    Vec_D     = (2.0:0.02:3).';
    Vec_R     = (2.2:0.02:3.3).';
    Vec_L0    = (2:0.02:2.8).';
    Vec_H0    = (40:2:60).';

    N_Phi   = numel(Vec_Phi);
    N_Theta = numel(Vec_Theta);
    N_D     = numel(Vec_D);
    N_R     = numel(Vec_R);
    N_L0    = numel(Vec_L0);
    N_H0    = numel(Vec_H0);

    RMS = zeros(N_Phi,N_Theta, N_D, N_R, N_L0, N_H0);
    prod(size(RMS))./1e6




    for I_D=1:1:N_D
        [I_D, N_D]
        for I_Phi=1:1:N_Phi
            [Lat_Axis,Lon_Axis] = reckon(Lat_NCP,Lon_NCP,Vec_D(I_D),Vec_Phi(I_Phi));
            for I_L0=1:1:N_L0
                for I_H0=1:1:N_H0
                    Lt = L - Vec_L0(I_L0);
                    Ht = H - Vec_H0(I_H0);
                    [Lat_Poin,Lon_Poin] = reckon(Lat_Axis,Lon_Axis,Lt,Ht);
                    for I_R=1:1:N_R
                        for I_Theta=1:1:N_Theta
                            R = Vec_R(I_R);
                            Theta = Vec_Theta(I_Theta);
                            [Lat_Tel,Lon_Tel]   = reckon(Lat_Poin,Lon_Poin,R,Theta);

                            Dist = celestial.coo.sphere_dist_fast(Lon_Tel./RAD,Lat_Tel./RAD,HA_Tel./RAD,Dec_Tel./RAD);
                            RMS(I_Phi,I_Theta, I_D,I_R, I_L0,I_H0) = std(Dist);
                        end
                    end
                end
            end
        end
    end



    size(RMS)                
    [M,I]=Util.stat.minnd(RMS)

    save I.mat M I 

    clear RMS


else


    VecPhi   = (0:90:330).';
    VecH0    = (0:10:350).';
    VecTheta = (0:10:350).';

    N_Phi   = numel(VecPhi);
    N_H0    = numel(VecH0);
    N_Theta = numel(VecTheta);


    MinFval=Inf;
    for I_Phi=1:1:N_Phi
        for I_H0=1:1:N_H0
            for I_Theta=1:1:N_Theta
                [X,Fval]=Util.fit.fminsearch_my({@LossFun,HA_Axis,Dec_Axis,HA_Tel,Dec_Tel},[2 VecPhi(I_Phi) 2 VecH0(I_H0) 2 VecTheta(I_Theta)]); 
                %[VecPhi(I_Phi), VecH0(I_H0), VecTheta(I_Theta), Fval]
                if Fval<MinFval
                    MinFval = Fval
                    Sol     = [VecPhi(I_Phi), VecH0(I_H0), VecTheta(I_Theta), Fval, I_Phi, I_H0, I_Theta];
                end

            end
        end
    end

    Sol
    MinFval
end

end

function RMS=LossFun(Par,HA_Axis,Dec_Axis,HA_Tel,Dec_Tel)
%%

RAD = 180./pi;

L = 90 - Dec_Axis;
H = HA_Axis;

Lat_NCP = 90;
Lon_NCP = 0;


D     = Par(1);
Phi   = Par(2);
L0    = Par(3);
H0    = Par(4);
R     = Par(5);
Theta = Par(6);

[Lat_Axis,Lon_Axis] = reckon(Lat_NCP,Lon_NCP, D, Phi);

Lt = L - L0;
Ht = H - H0;
[Lat_Poin,Lon_Poin] = reckon(Lat_Axis,Lon_Axis,Lt,Ht);
  
[Lat_Tel,Lon_Tel]   = reckon(Lat_Poin,Lon_Poin,R,Theta);

Dist = celestial.coo.sphere_dist_fast(Lon_Tel./RAD,Lat_Tel./RAD,HA_Tel./RAD,Dec_Tel./RAD);
RMS  = imUtil.background.rstd(Dist).*RAD;

end