function Out=conv1_vargauss(Data,SigmaVec)
% Convolve an array, along the first dimension, with a variable gaussian.
% Package: Util
% Description: Convolve an array, along the first dimension (columns),
%              with a Gaussian with variable width as a function of
%              position.
%              All the Gaussians are normalized to unity between -Inf and
%              Inf.
% Input  : - Data to convolve (vector or matrix).
%          - Vector of Gaussian width with which to convolve the data.
%            The length of this vector should be equal to the length of the
%            data dimesnion along to convolve.
%            Alternatively this can be a function handle that return the
%            Gaussian sigma to use as a function of poistion index.
% Output : - Convolved output.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Util.filter.conv1_vargauss(rand(1024,5),@(X) 0.5+X./10);
% Reliable: 2
%--------------------------------------------------------------------------

[Length,SpatLen] = size(Data);
X = (1:1:Length).';
% evaluate the SigmaVec
if (isa(SigmaVec,'function_handle'))
    % SigmaVec is a function handle
    SigmaVec = SigmaVec(X);
end

% convolve each SpatL dim
Out = zeros(Length,SpatLen);

for Ix=1:1:Length
    Y = exp(-(Ix-X).^2./(2.*SigmaVec(Ix).^2))./(SigmaVec(Ix).*sqrt(2.*pi));
    Out(Ix,:) = Out(Ix,:) + nansum(bsxfun(@times,Data,Y),1);
end
