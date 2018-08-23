function [FlagCom,Flag]=array_select(Matrix,Operator,varargin)
% Select lines in an array which columns satisfay some criteria.
% Package: Util.array
% Description: Given a matrix, select lines which their column fullfill
%              a specific criteria. For example, the values in the second
%              columns are in some range or follow a specific criterion.
% Input  : - A matrix.
%          - An operator, @all or @any, which to use in order to combine
%            all the crietria. For example, if @all then all the crietria
%            will combined using and and operator.
%            If empty, use @all.
%          * Arbitrary number of arguments containing cell arrays of
%            criteria information.
%            This can be:
%            {Col Min Max} - where Col is a column index, and this column
%                            should be >Min and <Max.
%            {Col @Fun, [Pars]} -In this case the criteria is
%                            Fun(Matrix(:,Col),Pars), which return either
%                            false or true.
% Output : - A vector of flags (per line) indicating if line was fullfiling
%            (and/or) the list of criteria.
%          - A matrix of flags (per line, and column per criteria),
%            indicating if line was fullfiling each one of the criteria.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [FlagCom,Flag]=array_select(rand(1000,5),[],{1 0.8 0.9},{2 @gt 0.5});
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
    Operator = [];
end

if (isempty(Operator))
    Operator = @all;
end

Nlines = size(Matrix,1);
Narg   = numel(varargin);
Flag   = true(Nlines,max(1,Narg)); 
for Iarg=1:1:Narg
    % options are: {Col Min Max} or {Col @Fun} or {Col @Fun Val,...}
    Col = varargin{Iarg}{1};
    if (isnumeric(varargin{Iarg}{2}))
        % assume {Col Min Max}
        Flag(:,Iarg) = Matrix(:,Col)>varargin{Iarg}{2} & Matrix(:,Col)<varargin{Iarg}{3};        
    else
        % assume {Col @Fun} or {Col @Fun Val,...}
        Fun          = varargin{Iarg}{2};
        if (numel(varargin{Iarg})>2)
            Pars = varargin{Iarg}(3:end);
        else
            Pars = {};
        end
        Flag(:,Iarg) = Fun(Matrix(:,Col),Pars{:});
    end
end

FlagCom = Operator(Flag,2);

    