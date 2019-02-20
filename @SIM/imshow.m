function HIMAGE=imshow(Sim)
% Show the Sim image use imshow
% Package: @SIM
% Description: Show the images in a SIM object using imshow.m.
% Input  : - A SIM object.
%
% Output : - the handle to the image object created by imshow.
% License: GNU general public license version 3
% Tested : Matlab R2017b
%     By : vancky                    Nov 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: imshow(S)
%          H=imshow(S(1));
% Reliable: 2
%--------------------------------------------------------------------------
if numel(Sim)>1
    warning('Input more than one image, only show the first one!');
    Sim=Sim(1);
end
if  isempty(Sim.Im)
    error('No image in Sim,please check');
end
if nargout>0
HIMAGE=imshow(Sim.Im,[Sim.background.BackIm,Sim.background.BackIm+3*Sim.background.ErrIm]);% 3 sigma
else
    imshow(Sim.Im,[Sim.background.BackIm,Sim.background.BackIm+3*Sim.background.ErrIm],'InitialMagnification', 'fit');
end

