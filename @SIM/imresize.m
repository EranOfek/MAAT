function Sim=imresize(Sim,Scale,Method,varargin)
% Resize (scale) images in a SIM object using imresize.
% Package: @SIM
% Description: Resize (scale) images in a SIM object using imresize.m.
% Input  : - A SIM object.
%          - Scaling factor by which to resize the image.
%          - Interpolation method. See imresize.m for options.
%            Default is 'lanczos2'.
%          * Additional parameters to pass to ufun2sim.m.
% Output : - SIM object with resized images.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: imresize(S,2);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3)
    Method = 'lanczos2';
end


DefV.FunAddPar          = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Sim = ufun2sim(Sim,@imresize,varargin{:},'FunAddPar1',Scale,'FunAddPar',{Method,InPar.FunAddPar{:}});
