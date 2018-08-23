function [Cheb]=chebyshev_polyrep(Kind)
% SHORT DESCRIPTION HERE
% Package: Util
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cheb]=Util.fun.chebyshev_polyrep(2)
% Reliable: 
%--------------------------------------------------------------------------

switch Kind
    case 1
        K = 0; 
        Cheb{K+1} = [1];   % 0
        K = K + 1;
        Cheb{K+1} = [1 0];  % 1
        K = K + 1;
        Cheb{K+1} = [2 0 -1];  % 2
        K = K + 1;
        Cheb{K+1} = [4 0 -3 0];  % 3
        K = K + 1;  
        Cheb{K+1} = [8 0 -8 0 1];  % 4
        K = K + 1;
        Cheb{K+1} = [16 0 -20 0 5 0];  % 5
        K = K + 1;
        Cheb{K+1} = [32 0 -48 0 18 -1];  % 6
        K = K + 1;
        Cheb{K+1} = [64 0 -112 0 56 0 -7 0];  % 7
    case 2
        K = 0; 
        Cheb{K+1} = [1];   % 0
        K = K + 1;
        Cheb{K+1} = [2 0];  % 1
        K = K + 1;
        Cheb{K+1} = [4 0 -1];  % 2
        K = K + 1;
        Cheb{K+1} = [8 0 -4 0];  % 3
        K = K + 1;  
        Cheb{K+1} = [16 0 -12 0 1];  % 4
        K = K + 1;
        Cheb{K+1} = [32 0 -32 0 6 0];  % 5
        K = K + 1;
        Cheb{K+1} = [64 0 -80 0 24 -1];  % 6
        K = K + 1;
        Cheb{K+1} = [128 0 -192 0 80 0 -8 0];  % 7
    otherwise
        error('Unknown Kind option');
end
