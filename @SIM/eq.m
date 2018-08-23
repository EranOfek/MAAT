function SimOut=eq(Sim1,Sim2,varargin)
% equal (==) operator on SIM objects.
% Package: @SIM
% Description: equal (==) operator on SIM objects.
% Input  : - A SIM object.
%          - A SIM object or a scalar.
%          * Additional arguments to pass to bfun2sim.m
% Output : - A SIM object in which the images contains the results
%            of the == operator.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim1==5
% Reliable: 2
%--------------------------------------------------------------------------


%SimOut = sim_imarith('In1',Sim1,'In2',Sim2,'Op','==',varargin{:});
SimOut=bfun2sim(Sim1,Sim2,@eq,varargin{:});
