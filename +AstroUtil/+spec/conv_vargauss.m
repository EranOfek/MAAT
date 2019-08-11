function [VG]=conv_vargauss(Spec,ResVec,varargin)
% Convolution with a Gaussian with a wavelngth dependent width
% Package: AstroUtil
% Description: Convolution with a Gaussian with a wavelngth dependent
%              width.
% Input  : - A spectrum [Wavelength, Intensity]
%          - A resolution vector, with the same length of the spectrum.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ResolutionType' - 'lambda' | 'R'. Default is 'R'.
%            'InterpMethod'   - Default is 'cubic'
%            'FilterSizeMaxSigma' - Default is 3.
% Output : - Convolved spectrum.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Conv]=AstroUtil.spec.conv_vargauss
% Reliable: 
%--------------------------------------------------------------------------



DefV.ResolutionType       = 'R';
DefV.InterpMethod         = 'cubic';
DefV.FilterSizeMaxSigma   = 3;
DefV.ColWave              = 1;
DefV.ColInt               = 2;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (nargin==0)
    % simulation mode
    Wave = (4000:1:9000)';
    Int  = zeros(size(Wave));
    
    Int(1000) = 100;
    Int(3000) = 100;
    Spec = [Wave, Int];
    ResVec    = 100.*(1 + (Wave - min(Wave))./(3.*range(Wave)));
    
end

Wave = Spec(:,InPar.ColWave);
Int  = Spec(:,InPar.ColInt);

switch lower(InPar.ResolutionType)
    case 'lambda'
        SigmaWave = ResVec;
    case 'r'
        SigmaWave = Wave./ResVec;
    otherwise
        error('unkwnon ResolutionType option');
end




% %SWN=SigmaWave./max(SigmaWave);
% SWN = flip(min(SigmaWave)/SigmaWave);
% 
% Factor = range(Wave)./sum(SWN);
% SigNew = max(SigmaWave)./Factor;
% NewRes = Factor.*SWN;
% 
% NewWave = Wave(1) + [0; cumsum(NewRes)];
% NewInt  = interp1(Wave,Int,NewWave,InPar.InterpMethod,'extrap');
% 
% 
% % Gaussian Filter
% FilterSize = ceil(max(SigmaWave).*InPar.FilterSizeMaxSigma);
% 
% X      = (-FilterSize:1:FilterSize)';
% Y      = exp(-X.^2./(2.*SigNew.^2))./sqrt(2.*pi.*SigNew.^2);
% 
% Conv = conv(NewInt,Y,'same')



        
SWN=SigmaWave./max(SigmaWave);
Factor = range(Wave)./sum(SWN);
SigNew = max(SigmaWave)./Factor;
NewRes = Factor.*SWN;
NewWave = Wave(1);
i=1;
while NewWave(i)<Wave(end)
    y = interp1(Wave,NewRes,NewWave(i),'linear','extrap');
    temp = NewWave(i) +y;
    NewWave = [NewWave temp];
    i=i+1;
end
NewRes = interp1(Wave,NewRes,NewWave,'linear','extrap');
NewInt = interp1(Wave,Int,NewWave,'linear','extrap');
% Gaussian Filter
FilterSize = ceil(max(SigmaWave).*InPar.FilterSizeMaxSigma);

X = (-FilterSize:1:FilterSize)';
Y = exp(-X.^2./(2.*SigNew.^2))./sqrt(2.*pi.*SigNew.^2);

Conv = conv(NewInt,Y,'same');
VG = [NewWave.', Conv.'];
 

    
    