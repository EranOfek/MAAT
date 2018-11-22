function [Chi2,Dof,AngRad,AngRadErr]=chi2_bb_photometry(Par,Phot,Ebv)
% Given photometric observations calculate \chi^2 for BB with a given T.
% Package: AstroUtil.spec
% Description: Given photometric observations calculate \chi^2 for
%              a black body with a given effective temperature and angular
%              radius.
% Input  : - [T(K),AngRad(radians)].
%          - A structure array with photometry. An element per band, with
%            the following fields:
%            .Family
%            .Band   - Band name.
%            .System - System name.
%            .Mag    - Magnitude.
%            .MagErr - Magnitude error.
%          - Eb-v [mag].
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Chi2,Dof,AngRad,AngRadErr]=AstroUtil.spec.chi2_bb_photometry([11000 1],tRes)
% Reliable: 
%--------------------------------------------------------------------------


% Phot contains:
% Phot().Band
% Phot().System
% Phot().Mag
% Phot().MagErr


    
if (min(size(Par))==1)
    % assume vector of T and AngRad = 1
    T  = Par(:);
    Fac    = 1;
    AngRad = Fac.*ones(size(T));
else
    Fac    = 1;
    T      = Par(:,1);
    AngRad = Par(:,2);
end

AngRad = AngRad./(constant.pc.*10);  % 1 cm from 10 pc

Nt        = numel(T);
Nphot     = numel(Phot);
ObsMag    = zeros(Nphot,1);
ObsMagErr = zeros(Nphot,1);
Mag       = zeros(Nphot,Nt);
for I=1:1:Nphot
    ObsMag(I)    = Phot(I).Mag;
    ObsMagErr(I) = Phot(I).MagErr;
    [Mag(I,:)] = AstroUtil.spec.blackbody_mag_c(T,Phot(I).Family,Phot(I).Band,Phot(I).System,1,1./AngRad./constant.pc,Ebv).';
end

% column per Temperature, raw per obs
Offset = mean(ObsMag - Mag,1);
ErrOffset = std(ObsMag - Mag)./sqrt(Nphot-2);

AngRad = AngRad.*10.^(-0.2.*Offset(:));
AngRadErr = AngRad.*ErrOffset(:);   % linear approximation
%AngRadErr = AngRad.*log(10).*((1./10.^((2.*Offset(:))./5).*ErrOffset(:).^2)./25).^(1./2)

Chi2  = nansum((ObsMag - Mag-Offset).^2./ObsMagErr.^2,1);
Dof   = sum(~isnan(ObsMag)) - 2;

    
