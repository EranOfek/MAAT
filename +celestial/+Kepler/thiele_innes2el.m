function El=thiele_innes2el(varargin)
% Thiele-Innes to orbital elements
% Package: celestial.Kepler
% Description: Calculate the orbital elements from the Thiele-Innes
%              orbital elements.
% Input  : - A or a structure with the A,B,C,F,G,H fields.
%          - F
%          - B
%          - G
%          - C (optional)
%          - H (optional)
% Output : - Structure containing the 6 Thiele-Innes elements
%            (fields are A, B, C, F, G, H).
% See also: thiele_innes2.m
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: TI=celestial.Kepler.thiele_innes(1,1,1,1);
%          El=celestial.Kepler.thiele_innes2el(TI);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1),
   TI = varargin{1};
elseif (nargin==4),
   TI.A = varargin{1};
   TI.F = varargin{2};
   TI.B = varargin{3};
   TI.G = varargin{4};
elseif (nargin==6),
   TI.A = varargin{1};
   TI.F = varargin{2};
   TI.B = varargin{3};
   TI.G = varargin{4};
   TI.C = varargin{5};
   TI.H = varargin{6};
else
   error('Illegal number of input arguments');
end

AG  =  TI.A + TI.G; % = a.*(1+cos(i)).*cos(omega+Omega)
AmG =  TI.A - TI.G; % = a.*(1-cos(i)).*cos(omega-Omega)
BmF =  TI.B - TI.F; % = a.*(1+cos(i)).*sin(omega+Omega)
mBF = -TI.B - TI.F; % = a.*(1-cos(i)).*sin(omega-Omega)

% (A+G).^2 + (B-F).^2 = a.^2.*(1+cos(i)).^2
% (A-G).^2 +(-B-F).^2 = a.^2.*(1-cos(i)).^2
% a = (sqrt((A+G).^2 + (B-F).^2) + sqrt((A-G).^2 +(-B-F).^2))./2
El.a = (sqrt(AG.^2 + BmF.^2) + sqrt(AmG.^2 +mBF.^2))./2;

% omega+Omega = atan2(A+G, B-F)
% omega-Omega = atan2(A-G, -B-F)
OpO = atan2(BmF, AG);
OmO = atan2(mBF, AmG);

CosI = AG./(El.a.*cos(OpO)) - 1;
El.I = acos(CosI);

El.omega = 0.5.*(OpO + OmO);
El.Omega = 0.5.*(OpO - OmO);

% the sign of the inclination
if (isfield(TI,'C') && isfield (TI,'H')),
   % C = a.*sin(omega).*sin(i)
   % H = a.*cos(omega).*sin(i)

   SinI = TI.C./(El.a.*sin(El.omega));

   El.I = atan2(CosI,SinI);
end





