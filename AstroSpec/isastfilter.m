function Ans=isastfilter(Obj)
%--------------------------------------------------------------------------
% isastfilter function                                           AstroSpec
% Description: Check if object is AstFilter class (astronomical filter).
% Input  : - Object.
% Output : - true or false.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: isastfilter(F);
% Reliable: 2
%--------------------------------------------------------------------------

Ans = isa(Obj,'AstFilter');
