function [RMS]=fit_extinction(Spec,Temp,VecEbv,VecR);
%------------------------------------------------------------------------------
% fit_extinction function                                            AstroSpec
% Description: Given a spectrum and a template, find the best extinction
%              curve that should be applied to the spectrum in order
%              to best fit the template.
% Input  : - Spectrum [Wavelength[Ang], Flux].
%          - Template [Wavelength[Ang], Flux].
%          - Vector of E_{B-V} to test.
%          - Vector of R_{V} to test.
% Output : - Matrix of fit rms as a function of E_{B-V} and R_{V}.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jul 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------
BinSize = 60;

Nr     = length(VecR);
Nebv   = length(VecEbv);



RMS = zeros(Nr,Nebv);
for Ir=1:1:Nr,
   for Iebv=1:1:Nebv,

      Aw = optical_extinction(VecEbv(Iebv),'B','V',Spec(:,1)./1e4,'C',VecR(Ir));
      UnextSpec = [Spec(:,1), Spec(:,2).*10.^(+0.4.*Aw)];
      Factor    = nanmean([Temp(:,2)./UnextSpec(:,2)]);

      %Bin   = binning([Spec(:,1),UnextSpec(:,2).*Factor],BinSize);
      %Noise = median(Bin(:,3).*sqrt(Bin(:,7)));

      %Chi2(Ir,Iebv) = nansum(((UnextSpec(:,2).*Factor - Temp(:,2))./Noise).^2);
      RMS(Ir,Iebv)  = nanstd(UnextSpec(:,2).*Factor - Temp(:,2));
   end
end

%Dof = length(Spec(:,1));
