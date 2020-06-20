function [EventRate,Err]=threshold_fluctuations(Th,varargin)
% Estimate the number of filtered local maxima above/below some threshold.
% Package: ImUtil
% Description: Estimate the two-sided false alarm rate for source detection
%              in astronomical images. The function generate random images,
%              filter the images with a Gaussian with a specified FWHM,
%              and find connected local maxima above or below the user
%              specified threshold.
% Input  : - Detection threshold. DEfault is 5.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FWHM'    - Gaussian PSF FWHM [pix]. Default is 2.5.
%            'Nsim'    - Number of simulated images. Default is 1000. 
%            'ImSize'  - Image Size [pix]. Default is 1024.
%            'PixScale'- Pixel scale [arcsec/pix].
% Output : - Event rate above below threshold per deg^2.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [EventRate,Err]=ImUtil.Im.threshold_fluctuations(5,'Nsim',100);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==0)
    Th = 5;
end

DefV.FWHM                 = 2.5;
DefV.Nsim                 = 1000;
DefV.ImSize               = 1024;
DefV.PixScale             = 1;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Sigma = InPar.FWHM./2.235;
F     = Kernel2.gauss(Sigma,Sigma,0);
Norm  = sqrt(sum(F(:).^2));


C1=0;
for Isim=1:1:InPar.Nsim
    R=randn(InPar.ImSize,InPar.ImSize);
    FR=filter2(F,R)./Norm;
    T1=(imregionalmax(abs(FR)).*abs(FR))>Th;
    C1=C1+sum(T1(:));

end
EventRate = C1./(InPar.Nsim.*(InPar.ImSize./(3600./InPar.PixScale)).^2);  % 2.92 per deg^2
Err       = EventRate./sqrt(C1);

