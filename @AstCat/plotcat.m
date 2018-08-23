function H=plotcat(AstC,ColX,ColY,varargin)
% Plot two columns in an AstCat object against each other.
% Package: @AstCat
% Description: Plot two columns in an AstCat object against each other.
% Input  : - An AstCat object.
%          - Column index, column name of an arithmatic extression on
%            column names that will represent the X-axis.
%          - Column index, column name of an arithmatic extression on
%            column names that will represent the Y-axis.
%          * Arbitrary number of arguments to pass to the plot.m function.
% Output : - Handle for the points plotted for each AstCat element.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: plotcat(S1,'FLUX_PSF','PSF_CHI2CR-PSF_CHI2','.');
% Reliable: 2
%--------------------------------------------------------------------------

Ncat = numel(AstC);
H    = zeros(size(AstC));
for Icat=1:1:Ncat,
    % Get X and Y from catalog
    X = col_arith(AstC(Icat),ColX,'mat');
    Y = col_arith(AstC(Icat),ColY,'mat');
    H(Icat) = plot(X,Y,varargin{:});
    hold on;
end
hold off;

if (~ischar(ColX)),
    ColX = colind2name(AstC(1),ColX);
end
if (~ischar(ColY)),
    ColY = colind2name(AstC(1),ColY);
end
Hx = xlabel(regexprep(ColX,'_','\_'));
set(Hx,'FontSize',18,'Interpreter','Latex');
Hy = ylabel(regexprep(ColY,'_','\_'));
set(Hy,'FontSize',18,'Interpreter','Latex');
