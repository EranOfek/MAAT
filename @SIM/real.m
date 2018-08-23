function SimOut=real(Sim,varargin)
% Calculate the real function for each element in a SIM object.
% Package: @SIM
% Description: Calculate the real function for each element in a SIM
%              object.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: real(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@real,varargin{:});
