function SimOut=log10(Sim,varargin)
% Base 10 logarithm of elements of a SIM object
% Package: @SIM
% Description: Calculate the log10 for each element in a SIM object.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: log10(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@log10,varargin{:});
