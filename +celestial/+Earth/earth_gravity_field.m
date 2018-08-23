function [V,W,EGM]=earth_gravity_field(R,Phi,Lambda,varargin)
% Calculate the Earth gravity field for a set of locations.
% Package: celestial.Earth
% Description: Calculate the Earth gravity field for a set of locations.
%              For both rotating and non rotating Earth.
%              Mean IRTF pole is assumed.
% Input  : - Column vector of radii from the Earth center at which to
%            calculate the Earth potential [cm].
%          - Column vector Geocentric latitude of at which to
%            calculate the Earth potential [rad].
%            Note that these are Geocentric rather than Geodetic
%            coordinates.
%          - Column vector Geocentric longitude of at which to
%            calculate the Earth potential [rad].
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Model' - Earth gravitational potential model:
%                      'DGM-1S'       - DGM-1s model (deg 250; default)
%                      'EIGEN-6C3stat'- EIGEN-6C3stat model (deg 1949)
%                      More models are available from:
%                      http://icgem.gfz-potsdam.de/ICGEM/
%            'EarthPeriod' - The Earth rotation period, used for
%                      calculating the potential on the moving surface of
%                      the Earth. Default is 86164.09890369732 s.
% Output : - The non-rotating Earth potential at the specified
%            locations [cgs units].
%          - The rotating Earth potential at the specified
%            locations [cgs units].
%          - A structure contains the Earth Gravity Model used.
% References: http://icgem.gfz-potsdam.de/ICGEM/
%             Seidelmann  (1992)
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [V,W,EGM]=celestial.Earth.earth_gravity_field(6731e5,1,1);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.Model  = 'dgm-1s';
DefV.EarthPeriod = 86164.09890369732; %[s] of UT1

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

OmegaE      = 2.*pi./InPar.EarthPeriod;

R      = R.';
Phi    = Phi.';
Theta  = pi./2 - Phi;
Theta  = Theta;
Lambda = Lambda.';

switch lower(InPar.Model)
    case 'dgm-1s'
        EGM = Util.IO.load2('EGM_DGM-1S.mat');
    case 'eigen-6c3stat'
        error('This model is not formatted correctly');
        EGM = Util.IO.load2('EGM_EIGEN-6C3stat.mat');
    otherwise
        error('Unknown Model option');
end


A  = EGM.Radius;
GM = EGM.GM_cgs; 
MaxL = EGM.Cat(end,EGM.Col.l);
Nl   = size(EGM.Cat,1);
Nt   = length(Theta);
MatP = zeros(Nl,Nt);
Ind  = 1;
for L=0:1:MaxL
    Ind = Ind + L;
    %[Ind, Ind+L]
    MatP(Ind:Ind+L,:) = legendre(L,cos(Theta),'sch');
end



Sum = sum( bsxfun(@power, (A./R), EGM.Cat(:,EGM.Col.l)).* ...
          ( bsxfun(@times,EGM.Cat(:,EGM.Col.C),cos(EGM.Cat(:,EGM.Col.m)*Lambda)) + ...
            bsxfun(@times,EGM.Cat(:,EGM.Col.S),sin(EGM.Cat(:,EGM.Col.m)*Lambda)) ).*MatP);

V = GM./R .*Sum;
if (nargout>1)
    X = R.*cos(Phi).*cos(Lambda);
    Y = R.*cos(Phi).*sin(Lambda);
    W = V + 0.5.*OmegaE.^2.*(X.^2 + Y.^2);        
end


