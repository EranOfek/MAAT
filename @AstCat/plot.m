function H=plot(AstC,varargin)
% Plot columns of AstCat object.
% Package: @AstCat
% Description: Plot columns of AstCat object.
% Input  : - An AstCat object.
%          - A cell array of two column names, or a vector of two column
%           indices to plot one against another.
%          * Additional parameters to pass to plot.m
% Output : - Matrix of handles, one per AstCat element.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: H=plot(AstC,{'X','Y'},'.');
%          H=plot(AstC,[1 2],'k.');
%          H=plot(AstC);
% Reliable: 2
%--------------------------------------------------------------------------


Narg = numel(varargin);
if (Narg==0),
    Col = [1 2];
elseif (Narg==1),
    Col = varargin{1};
    varargin = {};
else
    Col = varargin{1};
    varargin = varargin(2:end);
end

if (isempty(Col)),
    Col = [1 2];
end

ColInd = colname2ind(AstC,Col);

IsHold = ishold;
Ncat = numel(AstC);
H    = zeros(size(Ncat));
for Icat=1:1:Ncat,
    H(Icat) = plot(AstC(Icat).Cat(:,ColInd(1)),AstC(Icat).Cat(:,ColInd(2)),varargin{:});
    hold on;
end

if (IsHold),
    hold on;
else
    hold off;
end
