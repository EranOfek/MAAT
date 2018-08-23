function [ListEdge,ListCenter]=image_blocks(Size,Block,Buffer,Behav)
% Define boundries (ccdsec) for image subsections.
% Package: ImUtil.Im
% Description: Given image size, block size or number of blocks divide
%              define the bounderies of sub blocks.
% Input  : - Image size [X,Y]. Default is [2048 4096].
%          - A two element vector defines the block size. A positive number
%            is interpreted as a block size, while negaive number is
%            interpreted as the number of blocks.
%            Default is [-1 -1].
%          - Buffer size to add around blocks. Default is 0.
%          - A parameter indicating what to do if a whole number of blocks
%            can not be fitted into the image. Options are:
%            'simple' - 
% Output : - Four column matrix of list of blocks in image
%            [Xmin, Xmax, Ymin, Ymax].
%          - Four column matrix of list of block centers [Xcen, Ycen].
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [ListEdge,ListCenter]=ImUtil.Im.image_blocks([2048 4096],[256 256],10,'simple');
% Reliable: 2
%--------------------------------------------------------------------------

Def.Block  = [-1 -1];
Def.Buffer = 0;
Def.Behav  = 'simple';
if (nargin==1)
    Block  = Def.Block;
    Buffer = Def.Buffer;
    Behav  = Def.Behav;
elseif (nargin==2)
    Buffer = Def.Buffer;
    Behav  = Def.Behav;
elseif (nargin==3)
    Behav  = Def.Behav;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isempty(Block))
    Block = Def.Block;
end

BlockSize = zeros(1,2);
if (Block(1)<0)
    % Block size X is negative
    % this is the number of blocks
    BlockSize(1) = Size(1)./abs(Block(1));
else
    BlockSize(1) = Block(1);
end

if (Block(2)<0)
    % Block size Y is negative
    % this is the number of blocks
    BlockSize(2) = Size(2)./abs(Block(2));
else
    BlockSize(2) = Block(2);
end

switch lower(Behav)
    case 'simple'
        Xmin = (1:BlockSize(1):Size(1)-BlockSize(1).*0.5).';
        Xmax = (BlockSize(1):BlockSize(1):Size(1)).';
        if (numel(Xmin)>numel(Xmax))
            Xmax = [Xmax; Size(1)];
        end
        if (Xmax(end)<Size(1))
            Xmax(end) = Size(1);
        end
        
        
        Ymin = (1:BlockSize(2):Size(2)-BlockSize(2).*0.5).';
        Ymax = (BlockSize(2):BlockSize(2):Size(2)).';
        if (numel(Ymin)>numel(Ymax))
            Ymax = [Ymax; Size(2)];
        end
        if (Ymax(end)<Size(2))
            Ymax(end) = Size(2);
        end
        
        
    otherwise
        error('Unknwon Behav option');
end

[Xmin,Ymin] = meshgrid(Xmin,Ymin);
[Xmax,Ymax] = meshgrid(Xmax,Ymax);

ListEdge   = [Xmin(:), Xmax(:), Ymin(:), Ymax(:)];
ListEdge(:,1) = ListEdge(:,1) - Buffer;
ListEdge(:,2) = ListEdge(:,2) + Buffer;
ListEdge(:,3) = ListEdge(:,3) - Buffer;
ListEdge(:,4) = ListEdge(:,4) + Buffer;
ListEdge(ListEdge(:,1)<0,1) = 1;
ListEdge(ListEdge(:,3)<0,3) = 1;
ListEdge(ListEdge(:,2)>Size(1),2) = Size(1);
ListEdge(ListEdge(:,4)>Size(2),4) = Size(2);

ListCenter = [0.5.*(ListEdge(:,1) + ListEdge(:,2)), 0.5.*(ListEdge(:,3) + ListEdge(:,4))];



    