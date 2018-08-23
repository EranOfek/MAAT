function [Alpha,RMS,Best]=chi2_spectra(Spectra,Template,RedshiftVec)
%--------------------------------------------------------------------------
% chi2_spectra function                                          AstroSpec
% Description: Perform a chi^2 fitting between a spectrum and a template
%              spectrum as a function of redshift.
% Input  : - Spectrum [Wavelength, specific flux].
%          - Rest frame Template [Wavelength, specific flux].
%          - Redshift vector (default is: [0]).
% Output : - Vector of parameters in which to multiply Temnplate flux
%            in order to match it to spectrum, for each redshift.
%          - Vector of RMS for each redshift.
%          - Structure containing the best fit parameters.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                   Oct 2009
%    URL: http://weizmann.ac.il/home/eofek/matlab/
% Example: S=get_spectra('Gal_E',[0 0],[3.1 3.1],0.5);
%          T=get_spectra('Gal_E',[0 0],[3.1 3.1],0.0);
%          [Alpha,RMS,Best]=chi2_spectra(S,T,(0:0.01:1).');
%          plot((0:0.01:1).',RMS)
% Reliable: 2
%--------------------------------------------------------------------------

Nr = length(RedshiftVec);

SampRate = 3;
Method = 'linear';

BestRMS = Inf;
RMS     = zeros(Nr,1).*NaN;
Alpha   = zeros(Nr,1);
for Ir=1:1:Nr,
   RedTemplate = shift_spec(Template, RedshiftVec(Ir), 'w','wi');

   [S_Spec,S_Temp]=eq_sampling(Spectra, RedTemplate, SampRate, Method);
   %Npt = size(S_Spec,1);
   
   Alpha(Ir) = S_Spec(:,2)\S_Temp(:,2);

   Resid   = S_Spec(:,2) - S_Temp(:,2)./Alpha(Ir);
   RMS(Ir) = std(Resid);

   if (RMS(Ir)<BestRMS),
      Best.S_Spec = S_Spec;
      Best.S_Temp = S_Temp;
      Best.Resid  = Resid;
      Best.RMS    = RMS(Ir);
      Best.Alpha  = Alpha(Ir);
      Best.Redshift = RedshiftVec(Ir);
   end
end
