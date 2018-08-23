function SimOut=mtimes(Sim1,Sim2,varargin)
% Matrix multiplication between SIM arrays (*).
% Package: @SIM
% Description: Matrix multiplication between SIM arrays (*), or SIM array
%              and a scalar.
% Input  : - A SIM object.
%          - A SIM object or a scalar.
%          * Additional arguments to pass to bfun2sim.m.
% Output : - A SIM object with the result multiplication.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim1*Sim2
% Reliable: 2
%--------------------------------------------------------------------------


%SimOut = sim_imarith('In1',Sim1,'In2',Sim2,'Op','*',varargin{:});
SimOut=bfun2sim(Sim1,Sim2,@mtimes,varargin{:});
