function OutSpec=spec_convolve_resolution(Spectrum,R,Filter)
%--------------------------------------------------------------------------
% spec_convolve_resolution function                                 ImSpec
% Description: Given a spectrum and a constant resolution term
%              (dLambda/Lambda), convolve a filter which width equal
%              to the resolution, with the spectrum. This function can
%              be used to degrade a spectrum to a specific resolution.
% Input  : - A two column spectrum [wavelength, intesnity].
%          - Resolution. Default is 500.
%          - Filter {'gauss','mean','median'}. Default is 'mean'.
%            In Gauss the filter half size is equal to the Gaussian
%            sigma.
% Output : - Convolved spectrum.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Spec=spec_convolve_resolution(Tran.StdAirlessSpec,200,'gauss');
% Reliable: 2
%--------------------------------------------------------------------------


Def.R      = 500;
Def.Filter = 'mean';
if (nargin==1),
    R      = Def.R;
    Filter = Def.Filter;
elseif (nargin==2),
    Filter = Def.Filter;
elseif (nargin==3),
    
else
    
end

InPar.InterpMethod = 'linear';

Col.Wave = 1;
Col.Flux = 2;
% change sampling to logarithmic
MinWave = min(Spectrum(:,Col.Wave));
MaxWave = max(Spectrum(:,Col.Wave));
Npoint  = size(Spectrum,1);

LogWave = logspace(log10(MinWave),log10(MaxWave),Npoint.*5);
InterpFlux = interp1(Spectrum(:,Col.Wave),Spectrum(:,Col.Flux),LogWave,InPar.InterpMethod);
FilterHalfSize = ceil(LogWave(1)./diff(LogWave(1:2))./R .*0.5);


switch lower(Filter)
    case 'mean'
        SmoothFlux = runmean1(InterpFlux,FilterHalfSize);
    case 'median'
        SmoothFlux = medfilt1(InterpFlux,FilterHalfSize);
    case 'gauss'
        X = (-FilterHalfSize.*3:1:FilterHalfSize.*3).';  % 3 sigma range
        Y = fun_gauss([1, 0, FilterHalfSize],X);
        Y = Y./trapz(X,Y);
        SmoothFlux = conv(InterpFlux,Y.','same');
    otherwise
        error('Unknown Filter option');
end

% return to original wavelength
OutSpec = Spectrum;
OutSpec(:,Col.Flux) = interp1(LogWave,SmoothFlux,OutSpec(:,Col.Wave),InPar.InterpMethod);


        

