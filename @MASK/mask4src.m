function Flag=mask4src(Mask,CellInd,Operation)
% Return bit mask for specific pixels in images (sources).
% Package: @MASK
% Description: Return bit mask for specific pixels in images (sources).
% Input  : - A MASK object with a single element.
%          - Either a vector of indices or a cell array of vectors of
%            indices. For each vector of indices the and/or/xor over all
%            mask values corresponding to the indices will be calculated.
%          - Operation: @Util.array.bitor_array, @Util.array.bitand_array
%            Default is @Util.array.bitor_array.
% Output : - A column vector of bit masks. An element per source
%            (i.e., element in the cell array).
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: CI = ImUtil.Im.find_within_radius_cell(Size,X,Y,Radius,Circle)
%          Flag=mask4src(Mask,CI);
% Reliable: 2
%--------------------------------------------------------------------------

MaskField   = 'Mask';

if (nargin==2)
    Operation = @Util.array.bitor_array;
end

if (numel(Mask)>1)
    error('MASK element must have a single element');
end

if (~iscell(CellInd))
    CellInd = {CellInd};
end
Ncell = numel(CellInd);
Flag  = zeros(Ncell,1);
for Icell=1:1:Ncell
    if (~isempty(Mask.(MaskField)))
        Vec = Mask.(MaskField)(CellInd{Icell});
        Flag(Icell) = Operation(Vec(:),1);
    end
end