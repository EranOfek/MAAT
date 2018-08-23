function Kappa=kappa_from_alpha(AlphaX,AlphaY,PixelSize);
%-----------------------------------------------------------------------------
% kappa_from_alpha function                                             glens
% Description:   Calculate the surface density of a lens
%                               numerically from the deflection field.
% Input  : - Matrix of deflection in X (AlphaX) - the deflection is in pixels,
%            assuming d_ls/d_s=1.
%          - Matrix of deflection in Y (AlphaY) - the deflection is in pixels,
%            assuming d_ls/d_s=1.
%          - PixelSize in arcsec (per pixel), default is [0.05 0.05] (HST-ACS).
% Output : - The surface density of the lens in units of the critical
%            density for lensing.
% Tested : Matlab 7.0
%     By : Eran O. Ofek               June 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: Kappa=kappa_from_alpha(AlphaX.*DlsDs,AlphaY.*DlsDs);
%-----------------------------------------------------------------------------
RAD = 180./pi;
ARCSEC_DEG = 3600;

if (nargin==2),
   PixelSize = [0.05 0.05];
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (length(PixelSize)==1),
   PixelSize = [PixelSize, PixelSize];
end

% convert spacing to radians.
PixelSize = PixelSize./(RAD.*ARCSEC_DEG);

[Dxx,Dxy] = gradient(AlphaX.*PixelSize(1),PixelSize(1),PixelSize(2));
[Dyx,Dyy] = gradient(AlphaY.*PixelSize(2),PixelSize(1),PixelSize(2));

Kappa     = 0.5.*(Dxx + Dyy);



%%
if (1==0),
DlsDs=ad_dist([0.68 1.734])./ad_dist(1.734);
DlsDs=1;
Kappa=kappa_from_alpha(AlphaX.*DlsDs,AlphaY.*DlsDs);

SM=get_constant('SolM','cgs');
C=get_constant('c','cgs');
G=get_constant('G','cgs');
Pc=get_constant('pc','cgs');
D=ad_dist(0.68).*ad_dist([0.68 1.734]).*Pc./ad_dist(1.734);
Const=2.*D.*4.*pi.*G./(C.^2);
SurfDen=2.*Kappa./Const .* (ad_dist(0.68)./206000.*Pc.*0.05).^2./SM;
[MatX,MatY]=meshgrid([1:1:1500]);
X= 720;
Y= 744;
Rad=sqrt((MatX-X).^2+(MatY-Y).^2);
Radius=[100:10:600]';
for J=1:1:length(Radius),  I=find(Rad<=Radius(J)); Mass(J)=sum(SurfDen(I)); end
semilogy(Radius.*0.05.*ad_dist(0.68)./206000./1000,Mass)
end
