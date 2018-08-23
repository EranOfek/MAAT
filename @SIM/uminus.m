function SimOut=uminus(Sim,varargin)
% Multiply a SIM array by -1 (unary -1)
% Package: @SIM
% Description: Multiply a SIM array by -1 (unary -1)
% Input  : - SIM class object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - The SIM array multiplied by -1
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: -Sim
% Reliable: 2
%--------------------------------------------------------------------------

%SimOut = sim_imarith('In1',Sim,'In2',-1,'Op','.*',varargin{:});            
SimOut=ufun2sim(Sim,@uminus,varargin{:});

