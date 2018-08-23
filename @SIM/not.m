function SimOut=not(Sim,varargin)
% Not operator (~) on the elements in a SIM object.
% Package: @SIM
% Description: Not operator (~) on a SIM arrays using ufun2sim.m
% Input  : - SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM class with the resulted not operation.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ~Sim
% Reliable: 2
%--------------------------------------------------------------------------

%SimOut = sim_ufun(Sim,'Op',@not,varargin{:});
SimOut=ufun2sim(Sim,@not,varargin{:});