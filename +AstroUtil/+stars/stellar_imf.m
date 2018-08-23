function [DnDm,CumFracN,TotalFracN,CumFracM,TotalFracM]=stellar_imf(Mass,FunType)
% The stellar initial mass function
% Package: AstroUtil.stars
% Description: Return the stellar initial mass function in a given 
%              mass range.
% Input  : - Vector of mass [solar mass].
%          - Intial mass function type:
%            'Salpeter'  - Salpeter (1955)
%            'Scalo'     - Scalo (1986)
%            'SM1979'    - Scalo & Miller (1979)
%            Or alternativel, a matrix containing the IMF:
%            [M1 M2 Power_Law_Index; M1 M2 ...]
% Output : - Vector of dN/dM for each mass.
%          - Cumulative number of stars from 0 to Mass.
%          - Total fraction of number of stars within the range
%            min(Mass) to max(Mass), relative to the total mass in
%            the 0..Inf range.
%          - Cumulative mass of stars from 0 to Mass.
%          - Total fraction of stellar mass within the range
%            min(Mass) to max(Mass), relative to the total mass in
%            the 0..Inf range.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DnDm,CumFracN,TotalFracN,CumFracM,TotalFracM]=AstroUtil.stars.stellar_imf(5,'Salpeter');
% Reliable: 2
%--------------------------------------------------------------------------
ColM1 = 1;
ColM2 = 2;
ColX  = 3;

if (ischar(FunType)==1),
   switch FunType
    case 'Salpeter'
       IMF = [0.10 125 1.35];
    case 'Scalo'
       IMF = [0.10 0.18 -2.60;
   	      0.18 0.42 0.01;
   	      0.42 0.62 1.75;
              0.62 1.18 1.08;
              1.18 3.50 2.50;
              3.5  125  1.63];
    case 'SM1979'
       IMF = [0.10 1.00 0.25;
      	      1.00 2.00 1.00;
              2.00 10.0 1.30;
              10.0 125  2.30];
    otherwise
       error('Unknown FunType Option');
   end
else
   % User supplied IMF
   IMF = FunType;
end   

%--- Find normalizing factors needed to make the mass function continuus ---
N        = size(IMF,1);
Norm     = zeros(N,1);
Norm(1)  = 1;
for I=2:1:N,
   X_0  = IMF(I-1,ColX);
   X_1  = IMF(I,ColX);
   M1_0 = IMF(I-1,ColM1);
   M2_0 = IMF(I-1,ColM2);
   M1_1 = IMF(I,ColM1);
   M2_1 = IMF(I,ColM2);

   dNdM2 = Norm(I-1).*M2_0.^(-(1+X_0));

   dNdM1 = M1_1.^(-(1+X_1));

   Norm(I) = dNdM2./dNdM1;
end


%--- Integrate total mass within the normalized mass function ---
TotN = zeros(N,1);
TotM = zeros(N,1);

for I=1:1:N,
   X  = IMF(I,ColX);
   M1 = IMF(I,ColM1);
   M2 = IMF(I,ColM2);

   TotN(I) = (-Norm(I)./X).*(M2.^(-X) - M1.^(-X));

   TotM(I) = (Norm(I)./(1-X)).*(M2.^(1-X) - M1.^(1-X));

end
TotalNum  = sum(TotN);
TotalMass = sum(TotM);

DnDm     = zeros(size(Mass));
CumFracN = zeros(size(Mass));
CumFracM = zeros(size(Mass));
for I=1:1:N,
   X  = IMF(I,ColX);
   M1 = IMF(I,ColM1);
   M2 = IMF(I,ColM2);
      
   J = find(Mass>=M1 & Mass<M2);

   DnDm(J) = Norm(I).*Mass(J).^(-(1+X));

   CumFracN(J) = (-Norm(I)./X).*(Mass(J).^(-X) - M1.^(-X)) + sum(TotN(1:1:I-1));

   CumFracM(J) = (Norm(I)./(1-X)).*(Mass(J).^(1-X) - M1.^(1-X)) + sum(TotM(1:1:I-1));

end
% normalize
CumFracN = CumFracN./TotalNum;
CumFracM = CumFracM./TotalMass;

[~,MinI]  = min(Mass);
[~,MaxI]  = max(Mass);
TotalFracN  = CumFracN(MaxI) - CumFracN(MinI);
TotalFracM  = CumFracM(MaxI) - CumFracM(MinI);
