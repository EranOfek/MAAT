function plot_hr(B,R,Berr,Rerr,Dist,Sym);
%--------------------------------------------------------------------
% plot_hr function                                          plotting
% Description: Given a set of B and R calibrated magnitudes, plot 
%              a color magnitude diagram with overlayed main
%              sequence line assuming a given distance. 
% Input  : - B mag. vector
%          - R mag. vector
%          - B error mag. vector.
%          - R error mag. vector
%          - distance vector. (in pc).
%          - symbol, default is 'o'.
% Output : B-R vs. R graph with main sequence on distances
%          defined by Dist.
% Tested : Matlab 5.1
%     By : Eran O. Ofek           August 1998
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------
if (nargin<6),
   Sym = 'o';
end

errorxy([B-R,R,sqrt(Berr.^2 + Rerr.^2),Rerr],[1 2 3 4],Sym)
hold on;
xlabel('B-R');
ylabel('R');
axis('ij')

% plot main sequence

M_R  = zeros(4,1);
M_BR = zeros(4,1);

M_BR(1) = 0;
M_BR(2) = 2.8;
M_BR(3) = 3.3;
M_BR(4) = 4.5;

M_R(1) = 2.8.*M_BR(1) + 0.7;
M_R(2) = 2.8.*M_BR(2) + 0.7;
M_R(3) = 6.9.*M_BR(3) - 10.78;
M_R(4) = 3.0.*M_BR(4) + 2.1;

% distance modulus
for I=1:1:length(Dist),
   DistMod = 5.*log10(Dist(I)) - 5;

   % instrumental magnitude
   Ins_m_r = M_R + DistMod;

   if (I ~= 1),
      plot(M_BR,Ins_m_r,'g-');
      ZAMS_Title = strcat('ZAMS - ',num2str(Dist(I)),' pc');
      ZAMS = text(0.5.*(M_BR(1)+M_BR(2)),0.5.*(Ins_m_r(1)+Ins_m_r(2))-OffsetY,ZAMS_Title);
      set(ZAMS,'fontsize',5,'rotation',-25);

   else

      plot(M_BR,Ins_m_r,'.-');

      % plot spectral type
      Sp_A_BR = 0.0;
      Sp_F_BR = 0.5;
      Sp_G_BR = 1.0;
      Sp_K_BR = 1.5;
      Sp_M1_BR = 2.8;
      Sp_M4_BR = 3.2;

      OffsetX = 0.08;
      OffsetY = -0.4;

      FontSize = 12;
      Tick_A  = text(Sp_A_BR+OffsetX, 2.8.*Sp_A_BR + 0.7+OffsetY,'A0');
      set(Tick_A,'fontsize',FontSize);
      Tick_F  = text(Sp_F_BR+OffsetX, 2.8.*Sp_F_BR + 0.7+OffsetY,'F0');
      set(Tick_F,'fontsize',FontSize);
      Tick_G  = text(Sp_G_BR+OffsetX, 2.8.*Sp_G_BR + 0.7+OffsetY,'G0');
      set(Tick_G,'fontsize',FontSize);
      Tick_K  = text(Sp_K_BR+OffsetX, 2.8.*Sp_K_BR + 0.7+OffsetY,'K0');
      set(Tick_K,'fontsize',FontSize);
      Tick_M1 = text(Sp_M1_BR+OffsetX, 2.8.*Sp_M1_BR + 0.7+OffsetY,'M1');
      set(Tick_M1,'fontsize',FontSize);
      Tick_M4 = text(Sp_M4_BR+OffsetX, 3.0.*Sp_M4_BR+2.1+OffsetY,'M4');
      set(Tick_M4,'fontsize',FontSize);

      ZAMS_Title = strcat('ZAMS - ',num2str(Dist(I)),' pc');
      ZAMS = text(0.5.*(M_BR(1)+M_BR(2)),0.5.*(Ins_m_r(1)+Ins_m_r(2))-OffsetY,ZAMS_Title);
      set(ZAMS,'fontsize',5,'rotation',-25);
   end
end


% plot WD main sequence





hold off;



