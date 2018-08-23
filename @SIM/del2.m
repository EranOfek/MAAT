function Sim=del2(Sim,varargin)
% Discrete Laplacian for images in a SIM object.
% Package: @SIM
% Description: Calculate the discrete laplacian (del2) of images in a
%              SIM object.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object with the gradient in the X direction.
%          - A SIM object with the gradient in the Y direction.
% See also: SIM/gradient.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Sim]=del2(S);
% Reliable: 2
%--------------------------------------------------------------------------

Sim = ufun2sim(Sim,@del2,varargin{:});
