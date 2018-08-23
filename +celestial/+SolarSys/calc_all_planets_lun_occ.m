function [OccAll,Pl]=calc_all_planets_lun_occ
% Lunar occultations of planets

Coo(1,:) = [34.85 32 0]; % Tel-Aviv
Coo(2,:) = [34.763 30.596 0.900];  % Mizpe Ramon
Coo(3,:) = [34.933 29.550 0]; % Eilat
Coo(4,:) = [35 32.817 0]; % Haifa
Coo(5,:) = [35.167 31.783 0.800];  % Jeroslem
Coo(6,:) = [35.567 33.2 0.300];  % Kiriat Shmona
Nc = size(Coo,1);

VecYear = (2010:5:2045).';
Ny = length(VecYear);

for Ic=1:1:Nc,
   OccAll{Ic} = [];
   for Iy=1:1:Ny,
      [Occ{Iy},Pl]=planets_lunar_occultations([1 1 VecYear(Iy)],[1 1 VecYear(Iy)+5],'all',Coo(Ic,:),2);
      OccAll{Ic} = [OccAll{Ic}; Occ{Iy}]
   end
end

