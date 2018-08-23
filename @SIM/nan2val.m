function Sim=nan2val(Sim,Val,varargin)
% Replace NaN in SIM object images with a specified value.
% Package: @SIM
% Description: Replace NaN in SIM object images with a specified value.
% Input  : - A SIM object image array.
%          - A value that will replace all the NaNs in the images.
%            Default is 0. If empty then use default.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            to be passed to ufun2sim.m.
% Output : - A SIM object image array.
% License: GNU general public license version 3
% Tested : Matlab R2015b
% See also: nan2val.m
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=nan2val(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<2),
    Val = 0;
end

Sim = ufun2sim(Sim,@Util.array.nan2val,varargin{:},'FunAddPar',{Val});
