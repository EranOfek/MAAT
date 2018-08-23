function ImSize=imagesize(Sim,SizeFromImage,ImageField)
% Get the image size of a set of SIM images.
% Package: @SIM
% Description: Get the image size of a set of SIM images.
% Input  : - A SIM object.
%          - A flag indicating if to get size from image matrix size
%            (true), or from image header (false). Default is true.
%          - Image field name. Default is 'Im'.
% Output : - Two column matrix of image sizes [NAXIS1, NAXIS2].
%            I.e., [X, Y] size.
% See also: naxis.m
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ImSize=imagesize(Sim);
%          ImSize=imagesize(Sim,false);
% Reliable: 2
%--------------------------------------------------------------------------

Def.SizeFromImage = true;
Def.ImageField    = SIM.ImageField;

if (nargin==1)
    SizeFromImage = Def.SizeFromImage;
    ImageField    = Def.ImageField;
elseif (nargin==2)
    ImageField    = Def.ImageField;
elseif (nargin==3)
    % do nothing
else
    error('Illegal number of input arguments: imagesize(Sim,[SizeFromImage,ImageField])');
end



if (SizeFromImage)
    Nsim   = numel(Sim);
    ImSize = zeros(Nsim,2);
    for Isim=1:1:Nsim
        ImSize(Isim,:) = fliplr(size(Sim(Isim).(ImageField)));
    end
else
    ImSize = naxis(Sim);
end
    
