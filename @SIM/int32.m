function SimOut=int32(Sim,varargin)
% Convert all the image elements in a SIM array into int32.
% Package: @SIM
% Description: Convert all the image elements in a SIM array into int32.
%              By default this operates on the image field. In order to
%              run this on additional fields see ufun2sim.m
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object array with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: int32(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@int32,varargin{:});
