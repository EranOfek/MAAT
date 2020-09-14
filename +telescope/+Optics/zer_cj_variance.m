function [SigmaC,J,C]=zer_cj_variance(MaxJ,varargin)
% The expectency of std of the atmospheric Zernike coefficients.
% Package: telescope.Optics
% Description: Return the expectency of standard deviation of the
%              atmospheric Zernike c_j coefficients.
% Input  : - Maximum J. Default is 100.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'D'  - Telescope diamater [cm]. Default is 100.
%            'r0' - Fried length [cm]. Default is 10.
%            'Model' - The Zernike coefficient model to use:
%                      'zpl' - Power-law appropriate for Zenike functions
%                              (default).
%                      'pl'  - pure power law.
%            'Nrand' - Number of random relaizations of C to generate.
%                   Default is 100.
% Output : - Expectency of standard deviation of the atmospheric
%            Zernike c_j coefficients.
%          - Vector of J-s (i.e., (1:1:MaxJ)).
%          - Random realizations of C_j. Realization per raw.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AmpC,J,C]=telescope.Optics.zer_cj_variance(100);
% Reliable: 2
%--------------------------------------------------------------------------

Def.MaxJ = 100;
if (nargin==0)
    MaxJ = Def.MaxJ;
end
if (isempty(MaxJ))
    MaxJ = Def.MaxJ;
end


DefV.D       = 100;
DefV.r0      = 10;
DefV.Model   = 'zpl';
DefV.Nrand   = 100;
%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

J    = (1:1:MaxJ);
%AmpC = (InPar.D./InPar.r0).^(5./6).*J.^-(sqrt(3)./2);

switch lower(InPar.Model)
    case 'zpl'
        Ai=[1.03  0.582 0.134 0.111 0.088 0.065 0.059 0.053 0.046 0.04];
        Ac = 0.2944.*J.^(-sqrt(3)./2);
        Ac(1:1:min(length(Ai),MaxJ)) = Ai(1:1:min(length(Ai),MaxJ));
        
        SigmaC = [0, sqrt(diff(-Ac .*(InPar.D./InPar.r0).^(5./3)))];
    case 'pl'
        Ac = 0.2944.*J.^(-sqrt(3)./2);
        
        SigmaC = [0, sqrt(diff(-Ac .*(InPar.D./InPar.r0).^(5./3)))];
    
    case 'pl1'
        
        SigmaC = sqrt((InPar.D./InPar.r0).^(5./3).*J.^(-sqrt(3)./2));
        
    otherwise
        error('Unknown Model option');
end

        
if (nargout>2 && InPar.Nrand>0)
   C = randn(InPar.Nrand,MaxJ);
   C = bsxfun(@times,C,SigmaC);
end
