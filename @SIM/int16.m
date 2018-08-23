function SimOut=int16(Sim,varargin)
% Convert all the image elements in a SIM array into int16.
% Package: @SIM
% Description: Convert all the image elements in a SIM array into int16.
%              By default this operates on the image field. In order to
%              run this on additional fields see ufun2sim.m
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object array with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: int16(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@int16,varargin{:});
