function SimOut=acos(Sim,varargin)
% The inverse cosine for each element in a SIM object
% Package: @SIM
% Description: Calculate the inverse cosine for each element in a SIM
%              object.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: acos(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@acos,varargin{:});
