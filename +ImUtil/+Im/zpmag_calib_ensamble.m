function [Out,Stat]=zpmag_calib_ensamble(Mag,Err,JD,MagStar,ErrStar,varargin);
%---------------------------------------------------------------------------
% zpmag_calib_ensamble function                                      ImPhot
% Description: Relative photometric calibration using the ensamble method.
% Input  : - Matrix containing instrumental magnitudes for all stars in all
%            images, where each column represent a star and each row
%            represent an image (i.e., epoch). Example: Mag(Image,Star).
%            No NaN are allowed. 
%            If empty matrix (i.e., []) then generate random light
%            curve (test mode).
%          - Matrix containing instrumental magnitude errors for all stars,
%            where each column represent a star and each row represent
%            an image (i.e., epoch).
%          - Optional vector of JD.
%            If empty then assign running index to observations.
%          - Magnitude/s of star/s to measure in each image Mag(Image,Star).
%          - Error/s of star/s to measure in each image Mag(Image,Star).
% Output : - Structure containing:
%            .Mag  : Matrix of magnitudes brought to the same relative
%                    zero point.
%            .Err  : Matrix of errors for each star.
%          - Structure containing the best fit information:
% Tested : Matlab 7.6
%     By : Eran O. Ofek                     April 2009
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------

Def.JD     = [];

if (nargin==2),
   JD = Def.JD;
end

% Get number of stars and images:
[Nim,Nst]    = size(Mag);

if (isempty(JD)==1),
  JD = [1:1:Nim].';
end        

%--- convert all magnitudes to flux ---
ZP          = 25;
Flux        = 10.^(-0.4.*(Mag - ZP));
FluxErr     = -0.4.*log(10).*Flux.*Err;
FluxStar    = 10.^(-0.4.*(MagStar - ZP));
FluxErrStar = -0.4.*log(10).*FluxStar.*ErrStar;

ImageEnsambleFlux = sum(Flux,2);

Nstar     = size(MagStar,2);   % number of stars to measure
Nref      = size(Mag,2);       % number of reference stars
Nim       = size(Mag,1);       % number of images

% relative flux of stars to measure
Out.FluxStar = FluxStar./ImageEnsambleFlux(:,ones(1,Nstar));

Out.Flux = Flux./ImageEnsambleFlux(:,ones(1,Nref));
Out.Mag  = -2.5.*log10(Out.Flux);

Stat.Chi2 = NaN;
