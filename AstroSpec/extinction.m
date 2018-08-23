function A_3=extinction(Ebv,Family,Filter,R)
%--------------------------------------------------------------------------
% extinction function                                            AstroSpec
% Description: Given the E_{B-V} and the wavelength or filter name calculate
%              the extinction in magnitude.
%              The program works in the 0.1-2 micron range.
%              The program is using the Cardelli, Clayton, Mathis (1989)
%              model.
%              See also old version: optical_extinction.m.
% Input  : - E_{B-V} [mag]
%          - Filter effective wavelength [microns] or filter family name
%            (e.g., 'SDSS'). See get_filter.m for details.
%          - Filter name (e.g., 'r').
%          - R (scalar) Av/E(B-V), default value is 3.08.
% Output : - Extinction [mag]
% Tested : MATLAB 5.1
%     By : Eran O. Ofek                    Feb 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: optical_extinction.m
% Reference : Cardelli, Clayton, Mathis 1989 ApJ 345, 245
% Example: % to calculate The extinction in 'U' band given E(B-V)=0.1
%          A_U = extinction(0.1,'Johnson','U')
% Reliable: 1
%--------------------------------------------------------------------------
if (nargin==2),
   Filter = [];
   R     = 3.08;
elseif (nargin==3),
   R     = 3.08;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isnumeric(Family)),
   L = Family;
else
   F = get_astfilter(Family,Filter);
   L = [F.eff_wl]./10000; 
end


X = 1./L;

N = length(L);
if (length(R)==1),
   R = R.*ones(N,1);
end


A = zeros(size(X)).*NaN;
B = zeros(size(X)).*NaN;

% 0.3-1.1
I = find(X>=0.3 & X<1.1);
A(I) =  0.574.*X(I).^1.61;
B(I) = -0.527.*X(I).^1.61;

% 1.1-3.3
I = find(X>=1.1 & X<3.3);
Y = X(I) - 1.82;
A(I) = 1 + 0.17699.*Y - 0.50447.*Y.^2 - 0.02427.*Y.^3 + 0.72085.*Y.^4 + 0.01979.*Y.^5 - 0.77530.*Y.^6 + 0.32999.*Y.^7;
B(I) = 1.41338.*Y + 2.28305.*Y.^2 + 1.07233.*Y.^3 - 5.38434.*Y.^4 - 0.62251.*Y.^5 + 5.30260.*Y.^6 - 2.09002.*Y.^7;

% 5.9-8
I = find(X>=5.9 & X<8);
Fa = -0.04473.*(X(I)-5.9).^2 - 0.009779.*(X(I)-5.9).^3;
Fb = 0.2130.*(X(I)-5.9).^2 + 0.1207.*(X(I)-5.9).^3;
A(I)  =  1.752 - 0.316.*X(I) - 0.104./((X(I)-4.67).^2 + 0.341) + Fa;
B(I)  = -3.090 + 1.825.*X(I) + 1.206./((X(I)-4.62).^2 + 0.263) + Fb;

% 3.3-5.9
I = find(X>=3.3 & X<5.9);
A(I) =  1.752 - 0.316.*X(I) - 0.104./((X(I)-4.67).^2 + 0.341);
B(I) = -3.090 + 1.825.*X(I) + 1.206./((X(I)-4.62).^2 + 0.263);

% 8-10
I = find(X>=8 & X<=10);
A(I) = -1.073 - 0.628.*(X(I)-8) + 0.137.*(X(I)-8).^2 - 0.070.*(X(I)-8).^3;
B(I) = 13.670 + 4.257.*(X(I)-8) - 0.420.*(X(I)-8).^2 + 0.374.*(X(I)-8).^3;


A_L3 = A + B./R;



Const = A_L3./Ebv;
%sprintf('\n The relation constant is : A(3)=E(1-2) x %4.2f \n', Const)

%Factor = (2.512.^E - 1)./(2.512.^E_L12 - 1);
%A_3    = 2.5.*log10(((2.512.^A_L3 - 1).*Factor) + 1);
A_3    = Ebv.*Const;
A_3    = A_L3.*R.*Ebv;

