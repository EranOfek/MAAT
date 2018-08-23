function [SimOut]=atan2(Sim1,Sim2,varargin)
% Atan2 for elements in a SIM object
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
% Example: atan2(Sim1,Sim2);
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=bfun2sim(Sim1,Sim2,@atan2,varargin{:});