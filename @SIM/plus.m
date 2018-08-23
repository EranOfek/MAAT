function [SimOut]=plus(Sim1,Sim2,varargin)
% Add SIM arrays using bfun2sim
% Package: @SIM
% Description: Add SIM arrays using bfun2sim.m
% Input  : - First SIM object.
%          - Second SIM object, or a scalar.
%          * Additional arguments to pass to bfun2sim.m
% Output : - A SIM class with the resulted addition.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim1+Sim2
% Reliable: 2
%--------------------------------------------------------------------------

%SimOut = sim_imarith('In1',Sim1,'In2',Sim2,'Op','+',varargin{:});
SimOut=bfun2sim(Sim1,Sim2,@plus,varargin{:});