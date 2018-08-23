function Flag=nan_search(Sim,Fields)
% Search for SIM images that contains NaNs
% Package: @SIM
% Description: Go over all elements of a SIM object and look for NaNs.
% Input  : - A Sim object.
%          - Field to check. Default is
%          {'Im','BackIm','ErrIm','WeightIm','Mask','PSF'}.
% Output : - A matrix of size [number of images; number of fields] with
%            logical flag indicating if a NaN was found in each element.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Flag=nan_search(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<2)
    Fields = {'Im','BackIm','ErrIm','WeightIm','Mask','PSF'};
end
if (~iscell(Fields))
    Fields = {Fields};
end

Nf   = numel(Fields);
Nsim = numel(Sim);
Flag = zeros(Nsim,Nf);
for Isim=1:1:Nsim
    for If=1:1:Nf
        Flag(Isim,If) = any(isnan(Sim(Isim).(Fields{If})(:)));
    end
end
