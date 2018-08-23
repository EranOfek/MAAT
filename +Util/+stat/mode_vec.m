function [Mode,WMode,Step]=mode_vec(Vec,Npt)
%--------------------------------------------------------------------------
% mode_vec function                                              AstroStat
% Description: Calculate the mode of a vector using histogram.
%              The histogram step size chosen to contain some mean number
%              of points per bin.
%              The function uses only data within the 25% and 75%
%              percentile.
% Input  : - Vector, or an array that will be treated as a vector.
%          - Optional mean number of points per bin. Default is 200.
% Output : - Mode.
%          - The weighted mode. This is calculated by weighting the
%            histogram within the 25% and 75% percentile by the number of
%            counts in each bin.
%          - Step size of the histogram. This can be regarded as an
%            approximate resolutin of the mode.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mode,WMode]=mode_vec(randn(100000,1))
% Reliable: 2
%--------------------------------------------------------------------------

Def.Npt = 200;
Plow  = 0.25;
Phigh = 0.75;
if (nargin==1),
    Npt = Def.Npt;
end

Vec = Vec(~isnan(Vec(:)));
N   = numel(Vec);

Low  = quantile(Vec,Plow);
High = quantile(Vec,Phigh);

Nbin = ceil(N./Npt);

Step = (High-Low)./Nbin;

Edges    = (Low:Step:High);
[Ned]    = histcounts(Vec,Edges);
[~,MaxI] = max(Ned);
Mode     = Edges(MaxI);

if (nargout>1),
    % Calculate the weighted mode
    % assuming flat expectency
    BinCenter  = Edges(1:end-1)+Step.*0.5;
    WMode      = sum(BinCenter.*Ned)./sum(Ned);
end