function CC=ccf_fft2(Im1,Im2,SubMean,DivNorm)
% Cross correlation of 2D arrays usinf fft. 
% Package: Util.filter
% Description: Cross correlation of 2D arrays usinf fft. 
% Input  : - The first 2D array.
%          - The second 2D array.
%          - Subtract mean from each array {true|false}, default is false.
%          - Normalize the result {true|false}, default is false.
% Output : - Cross correlation of the two matrices.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Apr 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A=rand(1024,1024); CC=Util.filter.ccf_fft2(A,A);
% Reliable: 
%--------------------------------------------------------------------------

Def.SubMean = false;
Def.DivNorm = false;

if (nargin==2)
    SubMean = Def.SubMean;
    DivNorm = Def.DivNorm;
elseif (nargin==3)
    DivNorm = Def.DivNorm;
else
    % do nothing
end
    
   

[P1,Q1] = size(Im1);
[P2,Q2] = size(Im2);
P = max(P1,P2);
Q = max(Q1,Q2);

N = numel(Im1);

if (SubMean)
    M1 = mean(Im1(:));
    M2 = mean(Im2(:));
else
    M1 = 0;
    M2 = 0;
end

if (DivNorm)
    Norm = N.*std(Im1(:)).*std(Im2(:));
else
    Norm = 1;
end

CC = (ifft2(fft2(Im1-M1,P,Q).*conj(fft2(Im2-M2,P,Q)))) .* (1./Norm);
