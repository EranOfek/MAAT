function Var=load2(MatFile,varargin)
% Load a mat file into a variable
% Package: Util.IO
% Description: load a mat file containing a single variable to a variable
%              name (rather than a structure, like load.m).
%              If multiple variables are returned then will behave like
%              load.m
% Input  : - Mat file name.
%          * Additional parameters to pass to the load.m function.
% Output : - Variable name.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Jan 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------


Tmp = load(MatFile,varargin{:});
if (isstruct(Tmp))
   FN  = fieldnames(Tmp);
   if (length(FN)==1)
       Var = Tmp.(FN{1});
   else
       Var = Tmp;
   end
else
    Var = Tmp;
end

