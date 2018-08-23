function Sim=imrotate(Sim,Angle,varargin)
% Rotate the images in a SIM object using imrotate.
% Package: @SIM
% Description: Rotate the images in a SIM object using imrotate.m.
%              Note that the input images should be square, otherwise
%              output image size is not optimal.
% Input  : - A SIM object.
%          - Rotation angle in degrees in a counterclockwise direction
%            around its center point. This can be a scalar to be applied to
%            all aimages or an array with the same size as the SIM object.
%            In this case each element in the array will be used as the
%            rotation angle of the corresponding SIM element.
%          * Additional arguments to pass to ufun2sim.m.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MethodInterp' - Interpolation method. See imrotate.m for
%                             details. Default is 'blinear'.
%            'BBox'         - Bounding box to specify the size of the
%                             output image. See imrotate.m for
%                             details. Options are {'loose'|'crop'}.
%                             Default is 'loose'.
% Output : - Rotated image.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=imrotate(S1,45);
%          Sim=imrotate(S1(1:2),[45 47],'MethodInterp','linear','BBox','crop');
% Reliable: 2
%--------------------------------------------------------------------------

DefV.MethodInterp       = 'bilinear';
DefV.BBox               = 'loose';
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Sim = ufun2sim(Sim,@imrotate,varargin{:},'FunAddPar1',Angle,'FunAddPar',{InPar.MethodInterp,InPar.BBox});
