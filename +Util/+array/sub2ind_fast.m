function Ind=sub2ind_fast(Size,Y,X)
% sub2ind fast version for 2D matrices
% Description: A fast version of sub2ind for 2D arrays.
% Input  : - Array size [Y,X].
%          - Y positions (whole pixels)
%          - X positions (whole pixels)
% Output : - Linear index of position in array.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Sep 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Ind=sub2ind_fast([3 3],2,2)
% Reliable: 2
%--------------------------------------------------------------------------

Ind = uint32(Y + (X-1).*Size(1));