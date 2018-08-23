function NegCorr=sedm_wavelength_xcorr(K,Spec,Template,K1,FunDispLam,LowerK,UpperK)
%--------------------------------------------------------------------------
% sedm_wavelength_xcorr function                                      SEDM
% Description: Given two spectra and a distortion transformation
%              describing the relation between the wavelength of the
%              two spectra, return the negative of the cross-correlation
%              of the two spectra after applying the transformation.
% Input  : * FULL DESCRIPTION IS NOT AVAILABLE
%            This function is used by sedm_abswavecalib.m
% Output : - Negative of correlation.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: sedm_abswavecalib('SI',ArcSegmentsInfo(1130))
% Reliable: 2
%--------------------------------------------------------------------------



InPar.AngRes = 2;
InPar.Interp1Method = 'linear';
InPar.PolyDeg = 2;

%FunDispLam = @(K,Pix,K1) K1.*K(1).^(Pix-K(3)) + K(2).*(Pix-K(3)).^2; % + K(3);
W = FunDispLam(K,Spec(:,1),K1);
InterpSpec = interp1(W,Spec(:,2),Template(:,1),InPar.Interp1Method);

% subtract continuum from both spec and template
%ParIS = polyfit((1:1:length(InterpSpec)).',InterpSpec,InPar.PolyDeg);
%ParTE = polyfit(Template(:,1),Template(:,2),InPar.PolyDeg);
%InterpSpec = InterpSpec - polyval(ParIS,(1:1:length(InterpSpec)).');
%Template(:,2) = Template(:,2) - polyval(ParIS,Template(:,1));


FlagNN = ~isnan(InterpSpec) & ~isnan(Template(:,2));

CC = corrcoef(InterpSpec(FlagNN),Template(FlagNN,2));
if (size(CC,2)<2),
    NegCorr = 1; %return bad value
else
   NegCorr = -CC(1,2);
end

if (min(K-LowerK)<0 || max(K-UpperK)>0),
    NegCorr = 1;  % outof bound - return bad value
end

