function [AC]=sim2astcat(S)
% Copy a SIM object into AstCat object
% Package: @SIM
% Description: Copy a SIM object into AstCat object
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - An AstCat object
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AC]=sim2astcat(S);
% Reliable: 
%--------------------------------------------------------------------------


%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

N = numel(S);
AC = AstCat(size(S));

FN  = fieldnames(AC);
Nfn = numel(FN);

for I=1:1:N
    for Ifn=1:1:Nfn
        AC(I).(FN{Ifn}) = S(I).(FN{Ifn});
    end
end

    