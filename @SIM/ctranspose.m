function [obj]=ctranspose(Sim,varargin)
% ctranspose operator (') on a SIM array.
% Package: @SIM
% Description: ctranspose operator (') on a SIM array.
% Input  : - SIM object.
%          * Additional arguments to pass to sim_flip.m
% Output : - A SIM class with the resulted complex transpose operation.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim'
% Reliable: 2
%--------------------------------------------------------------------------

obj = sim_flip(Sim,'Op',@ctranspose,varargin{:});