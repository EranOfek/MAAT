function Vec=list2vec(varargin)
% Concatenate all vectors to a single vector.
% Package: Util.array
% Description: Given an arbitrary number of arguments, each containing
%              a vector, concatenate all vectors to a single vector.
% Input  : * Arbirtary number of arguments, each containing a vector.
%            The vectors can be either row or column vectors.
% Output : - A column vector, of all the vectors.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Apr 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Vec=list2vec(Out(I).SrcEnergyList);
% Reliable: 2
%--------------------------------------------------------------------------

Nin = length(varargin);
Vec = [];
for Iin=1:1:Nin,
    if (isempty(varargin{Iin})),
        % do nothing
    else
        Size = size(varargin{Iin});
        if (Size(1)>=Size(2)),
            Vec = [Vec; varargin{Iin}];
        else
            Vec = [Vec; varargin{Iin}.'];
        end
    end
end
