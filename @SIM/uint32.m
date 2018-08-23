function SimOut=uint32(Sim,varargin)
% Convert elements of a SiM object to uint32
% Package: @SIM
% Description: Convert all the image elements in a SIM array into uint32.
%              By default this operates on the image field. In order to
%              run this on additional fields see ufun2sim.m
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object array with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: uint32(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@uint32,varargin{:});
