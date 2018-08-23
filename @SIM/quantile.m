function SimOut=quantile(Sim,P,varargin)
%--------------------------------------------------------------------------
% quantile function                                             class/@SIM
% Description: Calculate the quantile of the values in each image in a
%              SIM object.
% Input  : - A SIM object.
%          - Percentile in fraction.
%          * Additional arguments to pass to ufun2scalar.m
% Output : - A matrix with the results - element per SIM element.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: quantile(S,5)
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2scalar(Sim,@quantile,varargin{:},'FunAddPar',{P});
