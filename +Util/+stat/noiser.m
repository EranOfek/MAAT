function Y=noiser(Y,varargin)
%--------------------------------------------------------------------------
% noiser function                                                AstroStat
% Description: Add multiple noise components to a vector.
% Input  : - Vector of values to which to add noise.
%          * Arbitrary number of ...,key,val, input arguments.
%            The following keywords are available:
%            'op'   - Add or multiply the noise {'+'|'*'}. Default is '+'.
%                     This can be either a char or a cell array. If cell
%                     array then each element refer to a specific noise
%                     element.
%            'noise'- Type of random numbers noise: {'norm','poiss'}.
%                     Default is 'norm'.
%                     This can be either a char or a cell array. If cell
%                     array then each element refer to a specific noise
%                     element.
%            'sigma'- Vector or scalar of sigma of normal distribution.
%                     Default is 1.
%            'sigmaL'- If this parameter is specified then use different
%                     noise for the upper and lower sigma of the normal
%                     distribution. While 'sigma' corresponds to the upper
%                     sigma, 'sigmaL' corresponds to the lower sigma.
%                     Default is empty.
%            'lambda'-Scalar or vectot of Lambda (expectency) of the
%                     Poisson distribution. Default is to use the first
%                     input argument (i.e., Y).
%                     
% Output : - Noisy vector
% Tested : Matlab R2011b                     
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% Y=noiser(rand(1000,1).*100,'op',{'+','+'},'noise',{'poiss','norm'},'sigma',6);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.op      = '+';
DefV.noise   = 'norm';
DefV.sigma   = 1;
DefV.sigmaL  = [];
DefV.lambda  = Y;
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

if (isempty(InPar.sigmaL)),
    InPar.sigmaL = InPar.sigma;
end

if (ischar(InPar.noise)),
    InPar.noise = {InPar.noise};
end
if (ischar(InPar.op)),
    InPar.op = {InPar.op};
end
Nnoise = length(InPar.noise);   % number of noise components

for Inoise=1:1:Nnoise,
    
    switch lower(InPar.noise{Inoise})
        case 'norm'
            Noise = randn(size(Y));
            SU    = (sign(Noise)+1).*0.5;
            SL    = -(sign(Noise)-1).*0.5;
            Noise = Noise.*SU.*InPar.sigma + Noise.*SL.*InPar.sigmaL;
        case 'poiss'
            Noise = poissrnd(Y) - Y;
        otherwise
            error('Unknown noise option');
    end
    
    switch lower(InPar.op{Inoise})
        case '+'
            Y = Y + Noise;
        case '*'
            Y = Y.*Noise;
        otherwise
            error('Unknown op option');
    end
end



