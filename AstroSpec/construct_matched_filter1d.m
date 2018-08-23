function Kernel=construct_matched_filter1d(Type,Pars,Norm)
%--------------------------------------------------------------------------
% construct_matched_filter1d function                            AstroSpec
% Description: Construct a 1D kernel.
% Input  : - Kernel type:
%            'gauss'  - Gaussian kernel. Parameters are [Half_Size, Sigma].
%            'tophat' - top hat kernel. Parameters are [Half_Size].
%            'triangle' - triangle kernel. Parameters are [Half_size].
%            'twohat' - Two top hats with zeros in between
%                       Parameters are [Half_size Half_gap],
%                       where the distance between the edges of the two
%                       top hats is 2*Half_size and the size of the gap
%                       between the top hats is 2*Half_gap.
%          - Parameters of kernel.
%          - Normalization of kernel. Default is 1.
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Kernel=construct_matched_filter1d('triangle',8)
%          Kernel=construct_matched_filter1d('twohat',[8 4]);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3),
    Norm = 1;
end

HalfSize = Pars(1);
X        = (-HalfSize:1:HalfSize).';
switch lower(Type)
    case 'gauss'
        % Gaussian kernel:
        Sigma = Pars(2);
        Kernel = exp(-X.^2./(2.*Sigma.^2));
    case 'tophat'
        % top hat kerenl
        Kernel = ones(size(X));
    case 'triangle'
        % triangle kernel
        X1 = (1:1:HalfSize).';
        X2 = (-HalfSize:1:-1).';
        Kernel = [0 + X1./(HalfSize+1); 1; 0 - X2./(HalfSize+1)];
    case 'twohat'
        % two top hats
        X1 = (1:1:HalfSize).';
        Kernel = ones(size(X1));
        Kernel(1:1:Pars(2)) = 0;
        Kernel = [flipud(Kernel); 0; Kernel];
    otherwise
        error('Unknown kernel option');
end

Kernel = Kernel.*Norm./trapz(X,Kernel);
