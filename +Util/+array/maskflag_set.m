function Mask=maskflag_set(Mask,Type,varargin)
%--------------------------------------------------------------------------
% maskflag_set function                                            General
% Description: Given a matrix or vector of bit masks, set specific
%              bits of specific indices.
% Input  : - Matrix of vector of bit masks. If empty will create a
%            new matrix which size is the size of the first Flag matrix.
%          - Matrix type to create (if not exist). Default is 'uint16'.
%            If empty use default.
%          * Arbitrary number of pairs of: ...,BitNumber,FlagBit,...
%            BitNumber is the index of the bit to set (e.g., 1),
%            while FlagBit is a matrix or vector of logicals (true|false)
%            of the indices in which to set the speciic bits.
% Output : - Output mask matrix.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: maskflag_check.m
% Example: Flag1=zeros(5,5,'uint16'); Flag2=ones(5,5,'uint16');
%          Flag3=zeros(5,5,'uint16'); Flag3(1,[2:3])=true;
%          Mask=maskflag_set([],[],1,Flag1,2,Flag2,4,Flag3);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Mask = [];
Def.Type = 'uint16';

if (isempty(Type)),
    Type  = Def.Type;
end

if (isempty(Mask)),
    Mask = zeros(size(varargin{2}),'uint16');
end

Narg = length(varargin);
for Iarg=1:2:Narg-1,
    if (isempty(find(varargin{Iarg+1}))),
        % do nothing
    else
       Mask(varargin{Iarg+1}==1) = bitset(Mask(varargin{Iarg+1}==1),varargin{Iarg});
    end
end

    
    