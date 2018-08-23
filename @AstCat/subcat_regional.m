function [SubCat,ListEdge,ListCenter]=subcat_regional(Cat,ImSize,BlockSize,BufferSize,Cols)
% Break an AstCat catalog into sub catalogs each for specific region.
% Package: @AstCat
% Description: Given an AstCat object with a single element break the
%              catalog into multiple catalogs each cover a sub region
%              in X and Y coordinates.
% Input  : - AstCat object with a single element.
%          - Image size [X,Y]. These are the maximal X and Y coordinates
%            in the catalog.
%          - The requsted block size (without buffer) [X, Y].
%            If 'full' then return the input as is.
%            Default is [1024 1024].
%          - Buffer size. Default is 200.
%          - Cell array of column names or vector of column indices of
%            X and Y columns. Default is [1 2].
% Output : - AstCat object with multiple elements. Element per regional
%            catalog.
%          - Matrix of [Xmin, Xmax, Ymin, Ymax] of all regions.
%          - Two column matrix of regions centers [X,Y].
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [SubCat]=subcat_regional(Cat,[2048 4096],[1024 1024],200,{'XWIN_IMAGE','YWIN_IMAGE'});
% Reliable: 2
%--------------------------------------------------------------------------

if (numel(Cat)>1)
    error('AstCat must have a single element');
end

if (ischar(BlockSize))
    switch lower(BlockSize)
        case 'full'
            SubCat = Cat;
            ListEdge = [];
            ListCenter = [];
        otherwise
            error('Unkown BlockSize option');
    end
else
    % blocks
    
    if all(BlockSize>=ImSize)
        ListEdge   = [1 ImSize(1) 1 ImSize(2)];
        ListCenter = ImSize.*0.5;
    else
        [ListEdge,ListCenter]=ImUtil.Im.image_blocks(ImSize,BlockSize,BufferSize);
    end
    Nblock = size(ListEdge,1);

    ColInd = colname2ind(Cat,Cols);

    if (SIM.issim(Cat))
        SubCat = SIM(Nblock,1);
    else
        if (AstCat.isastcat(Cat))
            SubCat = AstCat(Nblock,1);
        else
            error('Unknown input catalog class');
        end
    end

    for Iblock=1:1:Nblock
        SubCat(Iblock)  = select_ccdsec(Cat,ListEdge(Iblock,:),ColInd);
    end

end

