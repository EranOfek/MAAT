function varargout=index_outofbound(varargin)
% Remove from vector of indices valu8es which are out of bounds.
% Package: Util.array
% Description: Given a vector of indices and the allowed size of the array
%              will remove from the vector of indices any indices which
%              are out of bound.
% Input  : * Arbitrary number of pairs of arguments in which the first
%            of each pair is the vector of indices, and the second is
%            the allowed length.
% Output : * The number of output is equal to the number of pairs of
%            input arguments. Each output argument contains the new
%            vector of indices.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: I=index_outofbound((-1:1:7),6);
%          [I,J]=index_outofbound((-1:1:7),6,(-1:2:8),4);
% Reliable: 2
%--------------------------------------------------------------------------


Narg = length(varargin);
varargout = cell(1,Narg./2);
Counter = 0;
for Iarg=1:2:Narg-1,
    Ind = varargin{Iarg};
    Len = varargin{Iarg+1};
    
    Counter  = Counter + 1;
    varargout{Counter} = Ind(Ind>0 & Ind<=Len);
end
    