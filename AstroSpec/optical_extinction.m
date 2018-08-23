function A_3=optical_extinction(E,L1,L2,L3,Model,R);
%------------------------------------------------------------------------------
% optical_extinction function                                        AstroSpec
% Description: Given the E_{\lambda_1-\lambda_2} (e.g., E_{B-V}) calculate
%              the extinction in magnitude A_{\lambda_3}.
%              The program works in the 0.1-2 micron range.
%              The program is using the Cardelli, Clayton, Mathis (1989)
%              or Allen models.
%              \lambad can be specified using filters names:
%              Johnson: 'U','B','V','R','I','J','H','K'
%              POSS   : 'O','E'
%              SDSS   : 'u','g','r','i','z'
% Input  : - E_{\lambda_1-\lambda_2}
%          - \lambda_1, in microns, or by filter name.
%          - \lambda_2, in microns, or by filter name.
%          - \lmabda_3, in microns, or by filter name.
%          - Extinction law model:
%            'A' - Allen
%            'C' - Cardelli et al. (default).
%          - R (scalar) Av/E(B-V), default value is 3.08.
%            used only in the Cardelli et al. model
% Output : - A_{\lambda_3}
% Tested : MATLAB 5.1
%     By : Eran O. Ofek             Feb 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reference : Allen
%           : Cardelli, Clayton, Mathis 1989 ApJ 345, 245
% Example: % to calculate The extinction in 'U' band given E(B-V)=0.1
%          A_U = optical_extinction(0.1,'B','V','U','C',3.1);
% Reliable: 1
%------------------------------------------------------------------------------
if (nargin==4),
   Model = 'C';
   R     = 3.08;
elseif (nargin==5),
   R = 3.08;
elseif (nargin==6),
   % do nothing
else
   error('Illegal number of input arguments');
end
 
if (isstr(L1)==1),
   L1 = filter_eff_wavelength(L1);
end
if (isstr(L2)==1),
   L2 = filter_eff_wavelength(L2);
end
if (isstr(L3)==1),
   L3 = filter_eff_wavelength(L3);
end



% the extinction in the tables are normalized to A_{V}=1  !!!
switch Model
 case 'A'
    A_L1 = a_lambda_allen(L1);
    A_L2 = a_lambda_allen(L2);
    A_L3 = a_lambda_allen(L3);
 case 'C'
    A_L1 = a_lambda_cardelii(L1,R);
    A_L2 = a_lambda_cardelii(L2,R);
    A_L3 = a_lambda_cardelii(L3,R);
 otherwise
    error('Unknown extinction model');
end

E_L12 = A_L1 - A_L2;

Const = A_L3./E_L12;
%sprintf('\n The relation constant is : A(3)=E(1-2) x %4.2f \n', Const)

%Factor = (2.512.^E - 1)./(2.512.^E_L12 - 1);
%A_3    = 2.5.*log10(((2.512.^A_L3 - 1).*Factor) + 1);
A_3    = E.*Const;

A_3    = A_L3.*R.*E;

function W=filter_eff_wavelength(FilterName);
% find effective wavelength (micron) to filter name
%                      Johnson: 'U','B','V','R','I','J','H','K'
%                      POSS   : 'O','E'
%                      SDSS   : 'u','g','r','i','z'

switch FilterName,
 case 'U'
    W = 0.367;
 case 'B'
    W = 0.436;
 case 'V'
    W = 0.545;
 case 'R'
    W = 0.638;
 case 'I'
    W = 0.797;
 case 'J'
    W = 1.220;
 case 'H'
    W = 1.630;
 case 'K'
    W = 2.190;
 case 'O'
    W = 0.420;
 case 'E'
    W = 0.640;
 case 'u'
    W = 0.3557;
 case 'g'
    W = 0.4825;
 case 'r'
    W = 0.6261;
 case 'i'
    W = 0.7672;
 case 'z'
    W = 0.9097;
 otherwise
    error('Unkonw filter name');
end





function AAv = a_lambda_allen(L);
% given the wavelength in micron
% calc A/A_{\lambda} assuming R~3
% assuming Av=1;
% Reference : Allen

Lambda = [2.0
          1.0
          0.9
          0.67
          0.553
          0.50
          0.44
          0.40
          0.365
          0.333
          0.285
          0.250
          0.222
          0.200
          0.167
          0.143
          0.125
          0.111
          0.100];

A_Lam = [0.11
         0.38
         0.46
         0.74
         1.00
         1.13
         1.32
         1.45
         1.58
         1.69
         1.97
         2.30
         2.9
         2.8
         2.7
         3.0
         3.3
         3.7
         4.2];


AAv = interp1(Lambda, A_Lam, L);




function AAv = a_lambda_cardelii(L,R);
% given the wavelength in micron and R=Av/E(B-V)
% calc A/A_{\lambda}
% assuming Av=1;
% Reference : Cardelii et al.
if (nargin==1),
   R = 3.08;   % default
end


X = 1./L;

N = length(L);
if (length(R)==1),
   R = R.*ones(N,1);
end

AAv = zeros(N,1);
for I=1:1:N,
   if (X(I)>=0.3 & X(I)<=1.1),
      % IR
      A =  0.574.*X(I).^1.61;
      B = -0.527.*X(I).^1.61;
   elseif (X(I)>1.1 & X(I)<=3.3),
      % NIR/optical
      Y = X(I) - 1.82;
      A = 1 + 0.17699.*Y - 0.50447.*Y.^2 - 0.02427.*Y.^3 + 0.72085.*Y.^4 + 0.01979.*Y.^5 - 0.77530.*Y.^6 + 0.32999.*Y.^7;
      B = 1.41338.*Y + 2.28305.*Y.^2 + 1.07233.*Y.^3 - 5.38434.*Y.^4 - 0.62251.*Y.^5 + 5.30260.*Y.^6 - 2.09002.*Y.^7;
   elseif (X(I)>3.3 & X(I)<=8),
      % UV
      if (X(I)>=5.9),
         Fa = -0.04473.*(X(I)-5.9).^2 - 0.009779.*(X(I)-5.9).^3;
         Fb = 0.2130.*(X(I)-5.9).^2 + 0.1207.*(X(I)-5.9).^3;
      else
         Fa = 0;
         Fb = 0;
      end
      A =  1.752 - 0.316.*X(I) - 0.104./((X(I)-4.67).^2 + 0.341) + Fa;
      B = -3.090 + 1.825.*X(I) + 1.206./((X(I)-4.62).^2 + 0.263) + Fb;
   elseif (X(I)>8 & X(I)<=10),
      A = -1.073 - 0.628.*(X(I)-8) + 0.137.*(X(I)-8).^2 - 0.070.*(X(I)-8).^3;
      B = 13.670 + 4.257.*(X(I)-8) - 0.420.*(X(I)-8).^2 + 0.374.*(X(I)-8).^3;
   else
      error('Illegal wavelength, should be in range 0.1-3.33 micron');
   end

   AAv(I) = A + B./R(I);
end
