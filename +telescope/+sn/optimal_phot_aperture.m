function Res=optimal_phot_aperture(varargin)
% Optimal aperture radius for aperture photometry of a Gaussian PSF.
% Package: telescope.sn
% Description: Given a Gaussian symmetric PSF and image background and
%              readout noise. Estimate the radius of the photometric
%              aperture that maximize the S/N for an aperture photometry,
%              for a given target S/N. Return also the signal of the
%              source.
%              Calculate under the assumption that the number of pixels
%              in the background is continuous and not discrete.
%              Also calculate the S/N for a PSF photometry.
% Input  : * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Sigma'    - Gaussian sigma. Default is 3.
%            'TargetSN' - Tahrget S/N. Default is 3.
%            'B'        - Background [e-]. Default is 100.
%            'RN'       - Readout noise [e-]. Default is 5.
%            'Min_r'    - Minimum allowed aperture radius. Default is 1.
%            'Guess_r'  - Initial guess for optimal aperture. Default is 2.
%            'Thresh'   - Threshold for convergence. Default is 1e-3.
% Output : - A structure with the following fields:
%            .r  - radius of optimal aperture photometry.
%            .S  - Signal in source [e-] that will produce the target S/N.
%            .SN - S/N at convergence.
% See also: telescope.sn.sn_psf_phot.m
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Sep 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=telescope.sn.optimal_phot_aperture;
% Reliable: 2
%--------------------------------------------------------------------------

DefV.Sigma     = 3;
DefV.TargetSN  = 5;
DefV.B         = 100;
DefV.RN        = 5;
DefV.Min_r     = 1;
DefV.Guess_r   = 2;
DefV.Thresh    = 1e-3;
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});


FunFrac = @(r,sigma)             1 - exp(-r.^2./(2.*sigma.^2));
FunSN   = @(S,F,r,B,RN)          S.*F./sqrt(S.*F + pi.*r.^2.*(B + RN.^2));
FunS0   = @(TargetSN,F,r,B,RN)   TargetSN.^2./(2.*F) + TargetSN.*sqrt(TargetSN.^2.*F.^2 + 4.*pi.*r.^2.*F.^2.*(B+RN.^2) )./(2.*F.^2);



r1 = InPar.Guess_r;

% minimize S by varying r
Diff = Inf;
while (Diff>InPar.Thresh),
    r2 = r1.*1.1;
    r3 = r1.*1.2;
    F1  = FunFrac(r1,InPar.Sigma);
    F2  = FunFrac(r2,InPar.Sigma);
    F3  = FunFrac(r3,InPar.Sigma);
    
    S1  = FunS0(InPar.TargetSN,F1,r1,InPar.B,InPar.RN);
    S2  = FunS0(InPar.TargetSN,F2,r2,InPar.B,InPar.RN);
    S3  = FunS0(InPar.TargetSN,F3,r3,InPar.B,InPar.RN);
    % fit a 2nd deg polynomial
    Par = polyfit([r1;r2;r3],[S1;S2;S3],2);
    r1L = r1;
    % calculate the minimum of polynomial
    r1  = -Par(2)./(2.*Par(1));
    r1  = max(r1,InPar.Min_r);
    
    Diff = abs(r1L./r1 - 1);
end


Res.r  = r1;
Res.S  = S1;
Res.SN = FunSN(S1,F1,r1,InPar.B,InPar.RN);

% S/N for PSF photometry
Res.SN_PSFphot = telescope.sn.sn_phot_psf(Res.S,InPar.Sigma,InPar.B,InPar.RN);
