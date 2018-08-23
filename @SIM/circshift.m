function SimOut=circshift(Sim,K,varargin)
%--------------------------------------------------------------------------
% circshift function                                            class/@SIM
% Description: Apply the circshift function to each element in a SIM
%              object.
% Input  : - A SIM object.
%          - A scalar indicating by how much the elements in the array will
%            be shifted circularly along the first dimension.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: circshift(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@circshift,varargin{:},'FunAddPar',{K});
