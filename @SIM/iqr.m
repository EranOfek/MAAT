function SimOut=iqr(Sim,varargin)
% Calculate the interquartile range over all dimensions of SIM images.
% Package: @SIM
% Description: Calculate the iqr - interquartile range over all dimensions
%              for each image in a SIM object. This is the range between
%              the 25 to 75 percentiles.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2scalar.m
% Output : - A matrix with the results - element per SIM element.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: iqr(S)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2scalar(Sim,@iqr,varargin{:});
