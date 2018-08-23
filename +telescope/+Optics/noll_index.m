function [varargout]=noll_index(varargin)
% Return the Zernike index from the Noll index.
% Package: telescope.Optics
% Description: Calculate the Zernike indices given the Noll index and
%              visa versa.
% Input  : * Either j (the Noll index) or two parameters n and m.
% Output : * Either n and m or j (the Noll index).
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [n,m]=noll_index([1;6;8;100]);
%          [j]=noll_index([0;5;22],[0;-3;20]);
% Reliable: 2
%--------------------------------------------------------------------------
NollIndexList = cats.optics.NollIndexList;
%load NollIndexList.mat   % j,n,m % stored in matlab/data/General/

if (length(varargin)==1)
    % convert j to [n,m]
    varargout{1} = NollIndexList(varargin{1},2);
    varargout{2} = NollIndexList(varargin{1},3);
else
    % convert [n,m] to j
    N = numel(varargin{1});
    varargout{1} = zeros(N,1);
    for I=1:1:N
       
        J = NollIndexList(:,2)==varargin{1}(I) & NollIndexList(:,3)==varargin{2}(I);
        varargout{1}(I) = NollIndexList(J,1);
    end
end
