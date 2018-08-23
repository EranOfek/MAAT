function SubSim=subim_regional(Sim,BlockSize,BufferSize,varargin)
% Break a SIM image into sub images, with optional overlap.
% Package: @SIM
% Description: Given a SIM object with a single element, break the image
%              into sub-images which optional overlap.
% Input  : - A single element SIM object.
%          - Block size in pixels [I,J] (e.g., [1024 1024]).
%            Alternatively, if 'full' then return the original image.
%          - The overlap buffer size between sub images.
%            Default is 0.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - The SIM fields on which to execute the operator.
%                          Default is {'Im'}.
%            'TrimPar'   - The sub images are copied from the original
%                          image using SIM/trim_image.m.
%                          This is a cell array of additional key,val
%                          argument pairs to pass to trim_image.m
% Output : - A SIM object with the sub images.
%            The Original header is not copied.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: SubSim=subim_regional(S(2),[1024 1024],200);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField     = 'Im';

if (nargin==2)
    BufferSize = 0;
end

DefV.ExecField          = {ImageField};
DefV.TrimPar            = {};
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (numel(Sim)>1)
    error('SIM must have a single element');
end

ImSize = imagesize(Sim);
Nf     = numel(InPar.ExecField);

if (ischar(BlockSize))
    switch lower(BlockSize)
        case 'full'
            SubSim = Sim;
            ListEdge = [];
            ListCenter = [];
        otherwise
            error('Unkown BlockSize option');
    end
else
    % blocks
    
    [ListEdge,ListCenter]=ImUtil.Im.image_blocks(ImSize,BlockSize,BufferSize);
    Nblock = size(ListEdge,1);

    SubSim = simdef(Nblock,1);

    for Iblock=1:1:Nblock
        SubSim(Iblock)  = trim_image(Sim,ListEdge(Iblock,:),InPar.TrimPar{:},'ExecField',InPar.ExecField);
    end

end

