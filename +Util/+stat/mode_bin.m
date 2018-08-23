function [Mode]=mode_bin(Data)
%--------------------------------------------------------------------------
% mode_bin function                                              AstroStat
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

Plow   = 0.25;
Phigh  = 0.75;
NinBin = 100;
Dim    = 1;



[Ndata,Ncol] = size(Data);
Low   = quantile(Data,Plow,Dim);
High  = quantile(Data,Phigh,Dim);

Mode = zeros(Ncol,1);
for Icol=1:1:Ncol,
    BinSize = (High(Icol)-Low(Icol))./ceil(Ndata./NinBin)

    Edges = (Low:BinSize:High);
    Nbin  = histcounts(Data,Edges);
    bar(Edges(1:1:end-1),Nbin)
    
    [~,IndMax] = max(Nbin);
    Mode(Icol) = Edges(IndMax) + BinSize.*0.5;
end
