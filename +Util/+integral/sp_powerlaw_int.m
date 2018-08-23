function Res=sp_powerlaw_int(R,K,W1,R0,W2,Rlim)
%--------------------------------------------------------------------------
% sp_powerlaw_int function                                         General
% Description: Calculate the value the spherical integral and the line
%              integral of a broken power law of the form:
%              rho = K R^(-W1)               R<=R0
%              rho = K R0^(W2-W1) R^(-W2)    R>R0.
% Input  : - Matrix (or scalar) of X (e.g., radii) in which to calculate
%            the properties of the spherical broken power-law.
%          - The constant, K, that defines the power-law.
%            E.g., rho=Kr^-w (e.g., K is the mass-loading parameter).
%          - Power-law index of the first (inner) power-law.
%          - Radius at which the power-law is breaking.
%          - Power-law index of the second (outer) power-law.
%          - Limits for limit-integration [R1 R2].
%            If provided than will also calculate the line integral
%            between R1 and R2. Default is empty matrix.
% Output : - Structure with the following fields:
%            .Val   - Value of the function (e.g., Density) at the
%                     radius R.
%                     If K is something like the mass-loading parameter
%                     in gr*cm^-1 then
%                     in order to get particle density divide by
%                     <mu_p>*m_p
%            .Int   - Spherical integral of the function between
%                     radius 0 and R (e.g., mass).
%            .Col   - Line integral of the function between
%                     radius R and Infinity (e.g., column density).
%                     If K is something like the mass-loading parameter
%                     in gr*cm^-1 then
%                     in order to get particle density divide by
%                     <mu_p>*m_p, and multiply by sigmaT./(<mu_p>*m_p)
%                     in order to get the optical depth.
%            .ColLim  Line integral between R1 and R2. If R1 and R2 are
%                     not provided than will return NaN.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Apr 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=sp_powerlaw_int(logspace(10,18,100),9.7e16,2,4.4e15,4.36)
%          Res=sp_powerlaw_int(1,1e6,3.5,40,2.5,[0.001 1000]);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==5),
   Rlim = [];
end

I1 = find(R<=R0);
I2 = find(R>R0);

% value of the function at radii R
Val1 = K               .* R.^(-W1);     % R<=R0
Val2 = K.* R0.^(W2-W1) .* R.^(-W2);     % R>R0

Res.Val     = Val1;
Res.Val(I2) = Val2(I2);

% Spherical integral of the function between radius 0 and R (e.g., mass).
if (W1>3),
    Int1 = Inf.*ones(size(R));
    Int2 = Inf.*ones(size(R));
else
   Int1 = 4.*pi.*K./(3-W1) .* R.^(3-W1);   % R<=R0
   Int2 = 4.*pi.*K.*( (R0.^(3-W1))./(3-W1) - (R0.^(3-W1))./(3-W2) + (R0.^(W2-W1).*R.^(3-W2))./(3-W2) ); % R>R0
end

Res.Int     = Int1;
Res.Int(I2) = Int2(I2);

% Line integral of the function between radius R and Infinity (e.g., column density).
if (W2<1),
    Col1    = Inf.*ones(size(R));
    Col2    = Inf.*ones(size(R));
else
    Col1 = K.* R0.^(1-W1) ./(1-W1)  - K.* R.^(1-W1) ./(1-W1) - K.*R0.^(W2-W1) .* R0.^(1-W2)./(1-W2); % R<=R0
    Col2 = -K.* R0.^(W2-W1) .* R.^(1-W2)./(1-W2);     % R>R0
end

Res.Col     = Col1;
Res.Col(I2) = Col2(I2);

if (isempty(Rlim)),
   Res.ColLim = NaN;
else
   R1 = Rlim(1);
   R2 = Rlim(2);

   Fun1 = @(W1,K,R1,R2) K./(1-W1) .*(R2.^(1-W1) - R1.^(1-W1));
   Fun2 = @(W1,W2,K,R0,R1,R2) K./(1-W2).*R0.^(W2-W1).*(R2.^(1-W2) - R1.^(1-W2));
   if (R1<=R0 & R2<=R0),
      Res.ColLim = Fun1(W1,K,R1,R2);
   elseif (R1<R0 & R2>=R0),
      Res.ColLim = Fun1(W1,K,R1,R2) + Fun2(W1,W2,K,R0,R1,R2);
   elseif (R1>=R0 & R2>=R0),
      Res.ColLim = Fun2(W1,W2,K,R0,R1,R2);
   else
      error('Rlim is invalid');
   end
end
