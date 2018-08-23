function SimOut=transpose(Sim,varargin)
%Trranspose operator (.') on a SIM array.
% Package: @SIM
% Description: transpose operator (.') on a SIM array using ufun2sim.m
% Input  : - SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM class with the resulted transpose operation.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim.'
% Reliable: 2
%--------------------------------------------------------------------------

%SimOut = sim_flip(Sim,'Op',@transpose,varargin{:});
SimOut=ufun2sim(Sim,@transpose,varargin{:});

