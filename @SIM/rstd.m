function SimOut=rstd(Sim,varargin)
% Robust std over all dimensions for each image in a SIM object.
% Package: @SIM
% Description: Calculate the robust std over all dimensions for each image
%              in a SIM object. The robust std is calculated from the 50%
%              percentile.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2scalar.m
% Output : - A matrix with the results - element per SIM element.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: rstd(S)
% Reliable: 2
%--------------------------------------------------------------------------

if (isempty(varargin))
    SimOut=ufun2scalar(Sim,@Util.stat.rstd);
else
    SimOut=ufun2scalar(Sim,@Util.stat.rstd,'FunAddPar',varargin(:));
end