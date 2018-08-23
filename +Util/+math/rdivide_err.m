function [Ans,ErrAns]=rdivide_err(A,ErrA,B,ErrB)
% The rdivide operator including error propogation. 
% Package: Util.math
% Description: The rdivide operator including error propogation. 
% Input  : - A
%          - ErrA
%          - B
%          - ErrB
% Output : - A./B
%          - Error on A.*B. I.e., sqrt( (ErrA./B).^2 + (A.*ErrB./(B.^2)).^2 );
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Ans,ErrAns]=Util.math.rdivide_err(1,0.1,1,0.1)
% Reliable: 1
%--------------------------------------------------------------------------


Ans = A./B;
ErrAns = sqrt( (ErrA./B).^2 + (A.*ErrB./(B.^2)).^2 );
