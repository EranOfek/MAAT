function [SimX,SimY]=gradient(Sim,varargin)
% Calculate the gradient of images in a SIM object.
% Package: @SIM
% Description: Calculate the gradient of images in a SIM object.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object with the gradient in the X direction.
%          - A SIM object with the gradient in the Y direction.
% See also: SIM/del2.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [SimX,SimY]=gradient(S);
% Reliable: 2
%--------------------------------------------------------------------------

[SimX,SimY] = ufun2sim(Sim,@gradient,varargin{:});
