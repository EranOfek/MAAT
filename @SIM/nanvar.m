function SimOut=nanvar(Sim,varargin)
% Calculate the variance, and ignoring NaNs, of SIM image.
% Package: @SIM
% Description: Calculate the variance, and ignoring NaNs, over all
%              dimensions for each image in a SIM object.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2scalar.m
% Output : - A matrix with the results - element per SIM element.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: nanvar(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2scalar(Sim,@nanvar,varargin{:});
