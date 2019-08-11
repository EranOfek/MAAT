function [Chi2]=chi2_microlensing(Par,LC,Ndr,Fun,varargin)
% chi^2 between microlensing observations and model 
% Package: AstroUtil.microlensing
% Description: Calculate a \chi^2 for between a microlensing event
%              light curve and a theoretical model light curve
%              given microlensing parameters.
% Input  : - Parameters [T0, Beta, V, Alpha, BaseMag, SourceSize].
%            Where SourceSize is optional source size in Einstein
%            radius units. If given, the use psfs_microlensing.m
%            else use ps_microlensing.m.
%          - Light curve [JD, Mag, Err].
%          - Optional: Number of elements on the radius of the finite source,
%            in which to divide it when integrating its luminosity.
%          - Function, L=f(R,...), that given matrix R (radii) in units
%            of the Einstein radius and additional parameters, return the
%            luminosity of the finite source. The luminosity returned
%            by this function should be normalized such that the
%            integrated luminosity from the finite source is equal to
%            the base magnitude.
%            Default is to use a flat luminosity circular finite source.
%          * Arbitrary number of additional parameters to be passed to
%            the function that calculate the luminosity as a function
%            of position of the finite source.
% Output : - Chi^2
% Tested : Matlab 7.0
%     By : Eran O. Ofek                   October 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%-------------------------------------------------------------------------

import AstroUtil.microlensing.*

Col.JD  = 1;
Col.Mag = 2;
Col.Err = 3;

if (length(Par)==5),
   [Mag,Flux,Mu,Mu1,Mu2,U]=ps_microlens(LC(:,Col.JD),Par(1),Par(2),Par(3),Par(4),Par(5));
elseif (length(Par)==6),
   if (nargin==3),
      [Mag,Flux,Mu,Mu1,Mu2,U]=psfs_microlens(LC(:,Col.JD),Par(1),Par(2),Par(3),Par(4),Par(5),Par(6),Ndr);
   elseif (nargin==4),
      [Mag,Flux,Mu,Mu1,Mu2,U]=psfs_microlens(LC(:,Col.JD),Par(1),Par(2),Par(3),Par(4),Par(5),Par(6),Ndr,Fun);
   else
      [Mag,Flux,Mu,Mu1,Mu2,U]=psfs_microlens(LC(:,Col.JD),Par(1),Par(2),Par(3),Par(4),Par(5),Par(6),Ndr,Fun,varargin{:});
   end
else
   error('Length of Par vector should be 5 or 6');
end

Chi2 = sum( ((Mag - LC(:,Col.Mag))./LC(:,Col.Err)).^2 );

