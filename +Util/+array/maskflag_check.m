function FlagOn=maskflag_check(Mask,BinString,BitOp)
%--------------------------------------------------------------------------
% maskflag_check function                                          General
% Description: Given a matrix or vector of bit masks and a list of bits
%              to test, return true for indices in the matrix in which
%              one of the bits specified in the list of bits is open.
% Input  : - Matrix or vector of bit masks (in decimal representation).
%          - Either a string containing a binary number (e.g., '0101'),
%            or a decimal number representing the binary number
%            (e.g., 5 is equivalent to '0101'), or a vector of bits
%            (e.g., [0 1 0 1]).
%          - Operator {'or'|'and'} between the bits. Default is 'and'.
% Output : - Matrix or vector of logicals of the same size as the first
%            input argument in which the value is true if there are common
%            open bits between the bit mask and the binary number.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Oct 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FlagOn=maskflag_check([0 3 4],'100')
%          FlagOn=maskflag_check([0 3 4],[1 0 0])
%          FlagOn=maskflag_check([0 3 4],4)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2),
    BitOp = 'and';
end

if (ischar(BinString)),
    BinDec = bin2dec(BinString);
else
    if (numel(BinString)>1),
        BinString = sprintf('%d',BinString);
        BinDec = bin2dec(BinString);
    else
       BinDec = BinString;
    end
end

switch lower(BitOp)
    case 'and'
        if (BinDec==0),
            FlagOn = Mask == BinDec;
        else
            FlagOn = bitand(Mask,BinDec)~=0;
        end
    case 'or'
        if (BinDec==0),
            FlagOn = Mask == BinDec;
        else
            
            FlagOn = bitor(Mask,BinDec)~=0;
        end
            
    otherwise
        error('Unknown BitOp option');
end
