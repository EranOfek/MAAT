function SimOut=gpuArray(Sim,varargin)
%--------------------------------------------------------------------------
% gpuArray function                                                class/@SIM
% Description: Convert all the image elements in a SIM array into 
%              a gpuArray class so operations on the elements can be
%              run on the GPU.
%              By default this operates on the image field. In order to
%              run this on additional fields see ufun2sim.m
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object array with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S1=gpuArray(S);
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@gpuArray,varargin{:});
