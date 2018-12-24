function [Res,CatUM,Cat]=search_cat(AstC,varargin)
%--------------------------------------------------------------------------
% search_cat function                                        class/@AstCat
% Description: search_cat.m for AstCat object input.
%              See search_cat.m for details.
% Input  : - An AstCat object with a single element.
%          * Additional parameters to pass to search_cat.m.
% Output : * [Res,CatUM,Cat] - See search_cat.m for details.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,CatUM,Cat]=search_cat(AstC,[1 1]);
% Reliable: 2
%--------------------------------------------------------------------------

if (numel(AstC)>1)
    error('search_cat works on a single element AstCat object');
end
Def.ColX   = 1;
Def.ColY   = 2;

Ix = find(strcmpi(varargin(1:2:end-1),'ColX'));
if (isempty(Ix))
    ColX = Def.ColX;
else
    ColX = varargin{Ix.*2};
end

Iy = find(strcmpi(varargin(1:2:end-1),'ColY'));
if (isempty(Iy))
    ColY = Def.ColY;
else
    ColY = varargin{Iy.*2};
end

ColX = colname2ind(AstC,ColX);
ColY = colname2ind(AstC,ColY);

[Res,CatUM,Cat]=VO.search.search_cat(AstC.Cat(:,[ColX, ColY]),varargin{:});
