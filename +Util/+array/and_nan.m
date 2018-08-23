function B=and_nan(M1,M2);
% Logical function "and" for NaNs.
% Package: Util.array
% Description: Logical function "and" for NaNs. This function is similar
%              to "and" logical function, but NaNs are regarded as no
%              information using the following logical table:
%                    M1   M2    Result
%                    1    1     1
%                    1    0     0
%                    0    1     0
%                    0    0     0
%                    NaN  1     1
%                    NaN  0     0
%                    1    NaN   1
%                    0    NaN   0
%                    NaN  NaN   1
% Input  : - Matrix 1
%          - Matrix 2
% Output : - The result.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                  December 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

I    = find(isnan(M1)==1 | isnan(M2)==1); 
In   = find(isnan(M1)==1 & isnan(M2)==1); 
I1   = find(isnan(M1)==1); 
I2   = find(isnan(M2)==1); 
M1n  = M1;
M2n  = M2;
M1n(I1) = 0;
M2n(I2) = 0;
B       = M1n & M2n;
B(I)    = M1n(I) | M2n(I);   
B(In)   = 1;

